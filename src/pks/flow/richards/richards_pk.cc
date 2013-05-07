/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "Epetra_Import.h"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "upwinding.hh"
#include "Point.hh"
#include "matrix_mfd.cc"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content.hh"

#include "richards.hh"

#define DEBUG_RES_FLAG 0


namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards> Richards::reg_("richards flow");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
Richards::Richards(Teuchos::ParameterList& plist,
         const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist,solution),
    PKPhysicalBDFBase(plist, solution),
    coupled_to_surface_via_head_(false),
    coupled_to_surface_via_flux_(false),
    infiltrate_only_if_unfrozen_(false),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_wc_(false),
    modify_predictor_bc_flux_(false),
    upwind_from_prev_flux_(false),
    precon_wc_(false),
    niter_(0),
    dynamic_mesh_(false)
{
  // set a few parameters before setup
  plist_.set("primary variable key", "pressure");
  plist_.sublist("primary variable evaluator").set("manage communication", true);
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Richards::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupRichardsFlow_(S);
  SetupPhysicalEvaluators_(S);

  // debug cells
  if (coupled_to_surface_via_flux_ || coupled_to_surface_via_head_) {
    if (plist_.get<bool>("debug all surface cells", false)) {
      dc_.clear();
      Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh("surface");
      int ncells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
      for (int c=0; c!=ncells; ++c) {
        AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, c);
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);
        dc_.push_back(cells[0]);
      }
    }
  }
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void Richards::SetupRichardsFlow_(const Teuchos::Ptr<State>& S) {

  // Require fields and evaluators for those fields.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

#if DEBUG_RES_FLAG
  // -- residuals of various iterations for debugging
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->RequireField(namestream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
    S->RequireField(solnstream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  }
#endif

  // -- secondary variables, no evaluator used
  S->RequireField("darcy_flux_direction", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_velocity", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Create the absolute permeability tensor.
  int c_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (int c=0; c!=c_owned; ++c) {
    (*K_)[c].init(mesh_->space_dimension(),1);
  }

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();
  infiltrate_only_if_unfrozen_ = bc_plist.get<bool>("infiltrate only if unfrozen",false);
  bc_seepage_ = bc_factory.CreateSeepageFace();
  bc_seepage_->Compute(0.); // compute at t=0 to set up

  // how often to update the fluxes?
  std::string updatestring = plist_.get<std::string>("update flux mode", "iteration");
  if (updatestring == "iteration") {
    update_flux_ = UPDATE_FLUX_ITERATION;
  } else if (updatestring == "timestep") {
    update_flux_ = UPDATE_FLUX_TIMESTEP;
  } else if (updatestring == "vis") {
    update_flux_ = UPDATE_FLUX_VIS;
  } else if (updatestring == "never") {
    update_flux_ = UPDATE_FLUX_NEVER;
  } else {
    Errors::Message message(std::string("Unknown frequence for updating the overland flux: ")+updatestring);
    Exceptions::amanzi_throw(message);
  }

  // coupling
  // -- coupling done by a Neumann condition
  coupled_to_surface_via_flux_ = plist_.get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    S->RequireField("surface_subsurface_flux", name_)
        ->SetMesh(S->GetMesh("surface"))->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- coupling done by a Dirichlet condition
  coupled_to_surface_via_head_ = plist_.get<bool>("coupled to surface via head", false);
  if (coupled_to_surface_via_head_) {
    S->RequireField("surface_pressure");
    // override the flux update -- must happen every iteration
    update_flux_ = UPDATE_FLUX_ITERATION;
  }

  // -- Make sure coupling isn't flagged multiple ways.
  ASSERT(!(coupled_to_surface_via_flux_ && coupled_to_surface_via_head_));


  // Create the upwinding method
  S->RequireField("numerical_rel_perm", name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("numerical_rel_perm",name_)->set_io_vis(false);

  string method_name = plist_.get<string>("relative permeability method", "upwind with gravity");
  symmetric_ = false;
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
            "relative_permeability", "numerical_rel_perm", K_));
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            "relative_permeability", "numerical_rel_perm"));
    symmetric_ = true;
    Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwind_from_prev_flux_ = plist_.get<bool>("upwind flux from previous iteration", false);
    if (upwind_from_prev_flux_) {
      upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
              "relative_permeability", "numerical_rel_perm", "darcy_flux"));
    } else {
      upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
              "relative_permeability", "numerical_rel_perm", "darcy_flux_direction"));
    }
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            "relative_permeability", "numerical_rel_perm"));
    Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards FLow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, mesh_));
  matrix_->SetSymmetryProperty(symmetric_);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->InitPreconditioner();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  Teuchos::RCP<Operators::MatrixMFD> precon =
    Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, mesh_));
  set_preconditioner(precon);

  // wc preconditioner
  precon_wc_ = plist_.get<bool>("precondition using WC", false);

  // predictors for time integration
  modify_predictor_with_consistent_faces_ =
    plist_.get<bool>("modify predictor with consistent faces", false);
  modify_predictor_bc_flux_ = 
    plist_.get<bool>("modify predictor for flux BCs", false);
  modify_predictor_first_bc_flux_ = 
    plist_.get<bool>("modify predictor for initial flux BCs", false);
  modify_predictor_wc_ = 
    plist_.get<bool>("modify predictor via water content", false);
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Richards::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("water content evaluator");
  Teuchos::RCP<RichardsWaterContent> wc = Teuchos::rcp(new RichardsWaterContent(wc_plist));
  S->SetFieldEvaluator("water_content", wc);

  // -- Water retention evaluators
  // -- saturation
  Teuchos::ParameterList wrm_plist = plist_.sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);

  // -- rel perm
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";

  S->RequireField("relative_permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);
  wrms_ = wrm->get_WRMs();

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("molar_density_liquid");

  S->RequireField("viscosity_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("viscosity_liquid");

  // -- liquid mass density for the gravity fluxes
  S->RequireField("mass_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_density_liquid"); // simply picks up the molar density one.
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Richards::initialize(const Teuchos::Ptr<State>& S) {
  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  // debugggin cruft
#if DEBUG_RES_FLAG
  for (int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }
#endif

  // check whether this is a dynamic mesh problem
  if (S->HasField("vertex coordinate")) dynamic_mesh_ = true;

  // Initialize boundary conditions.
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S->GetFieldData("numerical_rel_perm",name_)->PutScalar(1.0);
  S->GetField("numerical_rel_perm",name_)->set_initialized();

  S->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
  S->GetField("darcy_flux", name_)->set_initialized();
  S->GetFieldData("darcy_flux_direction", name_)->PutScalar(0.0);
  S->GetField("darcy_flux_direction", name_)->set_initialized();
  S->GetFieldData("darcy_velocity", name_)->PutScalar(0.0);
  S->GetField("darcy_velocity", name_)->set_initialized();

  // initialize coupling terms
  if (coupled_to_surface_via_flux_) {
    S->GetFieldData("surface_subsurface_flux", name_)->PutScalar(0.);
    S->GetField("surface_subsurface_flux", name_)->set_initialized();
  }

  // absolute perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  matrix_->CreateMFDmassMatrices(K_.ptr());
  mfd_preconditioner_->CreateMFDmassMatrices(K_.ptr());
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Richards::commit_state(double dt, const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "Commiting state." << std::endl;

  niter_ = 0;

  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator(key_)->HasFieldChanged(S.ptr(), name_);
  update |= S->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S.ptr(), name_);

  if (update_flux_ == UPDATE_FLUX_TIMESTEP ||
      (update_flux_ == UPDATE_FLUX_ITERATION && update)) {

    Teuchos::RCP<const CompositeVector> rel_perm =
      S->GetFieldData("numerical_rel_perm");
    // update the stiffness matrix
    matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

    // derive fluxes
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
    Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", name_);
    matrix_->DeriveFlux(*pres, flux.ptr());
    AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
    flux->ScatterMasterToGhosted();
  }

  // As a diagnostic, calculate the mass balance error
#if DEBUG_FLAG
  if (S_next_ != Teuchos::null) {
    Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> wc0 = S_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux", name_);
    CompositeVector error(*wc1);

    for (int c=0; c!=error.size("cell"); ++c) {
      error("cell",c) = (*wc1)("cell",c) - (*wc0)("cell",c);

      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        error("cell",c) += (*darcy_flux)("face",faces[n]) * dirs[n] * dt;
      }
    }

    double einf(0.0);
    error.NormInf(&einf);

    // VerboseObject stuff.
    Teuchos::OSTab tab = getOSTab();
    *out_ << "Final Mass Balance Error: " << einf << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void Richards::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  if (update_flux_ == UPDATE_FLUX_VIS) {
    Teuchos::RCP<const CompositeVector> rel_perm =
      S->GetFieldData("numerical_rel_perm");
    // update the stiffness matrix
    matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

    // derive fluxes
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", name_);
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
    Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    matrix_->DeriveFlux(*pres, flux.ptr());
    AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
    flux->ScatterMasterToGhosted();
  }

  if (update_flux_ != UPDATE_FLUX_NEVER) {
    Teuchos::RCP<CompositeVector> darcy_velocity = S->GetFieldData("darcy_velocity", name_);
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("darcy_flux");
    matrix_->DeriveCellVelocity(*flux, darcy_velocity.ptr());

    S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& nliq_c = *S->GetFieldData("molar_density_liquid")
        ->ViewComponent("cell",false);

    Epetra_MultiVector& vel_c = *darcy_velocity->ViewComponent("cell",false);
    int ncells = vel_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      for (int n=0; n!=vel_c.NumVectors(); ++n) {
        vel_c[n][c] /= nliq_c[0][c];
      }
    }

  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Richards::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "  Updating permeability?";

  Teuchos::RCP<CompositeVector> uw_rel_perm = S->GetFieldData("numerical_rel_perm", name_);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  bool update_perm = S->GetFieldEvaluator("relative_permeability")
      ->HasFieldChanged(S, name_);

  // place n/mu on cells
  update_perm |= S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("viscosity_liquid")->HasFieldChanged(S, name_);
  const Epetra_MultiVector& n_liq = *S->GetFieldData("molar_density_liquid")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& visc = *S->GetFieldData("viscosity_liquid")
      ->ViewComponent("cell",false);

  // requirements due to the upwinding method
  if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    bool update_dir = S->GetFieldEvaluator("mass_density_liquid")
        ->HasFieldChanged(S, name_);
    update_dir |= S->GetFieldEvaluator(key_)->HasFieldChanged(S, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData("darcy_flux_direction", name_);

      // Create the stiffness matrix without a rel perm (just n/mu)
      for (int c=0; c!=uw_rel_perm->size("cell", false); ++c) {
        (*uw_rel_perm)("cell",c) = n_liq[0][c] / visc[0][c];
      }
      uw_rel_perm->ViewComponent("face",false)->PutScalar(1.);
      matrix_->CreateMFDstiffnessMatrices(uw_rel_perm.ptr());

      // Derive the pressure fluxes
      Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);
      matrix_->DeriveFlux(*pres, flux_dir.ptr());

      // Add in the gravity fluxes
      Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
      Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
      AddGravityFluxesToVector_(gvec.ptr(), uw_rel_perm.ptr(), rho.ptr(), flux_dir.ptr());
      flux_dir->ScatterMasterToGhosted();
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    // patch up the BCs -- move rel perm on boundary_faces into uw_rel_perm on faces
    const Epetra_Import& vandelay = mesh_->exterior_face_importer();
    const Epetra_MultiVector& rel_perm_bf =
        *rel_perm->ViewComponent("boundary_face",false);
    Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
    uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);

    // upwind
    upwinding_->Update(S);

    // patch up the surface, use 1
    //    if (coupled_to_surface_via_head_ || coupled_to_surface_via_flux_) {
    if (S->HasMesh("surface")) {
      Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh("surface");
      int ncells_surface = surface->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

      if (S->HasFieldEvaluator("unfrozen_fraction")) {
        S->GetFieldEvaluator("unfrozen_fraction")->HasFieldChanged(S.ptr(), name_);
        const Epetra_MultiVector& uf = *S->GetFieldData("unfrozen_fraction")
            ->ViewComponent("cell",false);
        for (int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face
          AmanziMesh::Entity_ID f =
              surface->entity_get_parent(AmanziMesh::CELL, c);

          // -- set that value to the unfrozen fraction to ensure we
          // -- don't advect ice
          uw_rel_perm_f[0][f] = uf[0][c];
        }
      } else {
        for (int c=0; c!=ncells_surface; ++c) {
          // -- get the surface cell's equivalent subsurface face
          AmanziMesh::Entity_ID f =
              surface->entity_get_parent(AmanziMesh::CELL, c);

          // -- set that value to 1
          uw_rel_perm_f[0][f] = 1.0;
        }
      }
    }
  }

  // Scale cells by n/mu
  if (update_perm) {
    for (int c=0; c!=uw_rel_perm->size("cell", false); ++c) {
      (*uw_rel_perm)("cell",c) *= n_liq[0][c] / visc[0][c];
    }
  }

  // debugging
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    *out_ << " " << update_perm << std::endl;
  }
  return update_perm;
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void Richards::UpdateBoundaryConditions_() {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "  Updating BCs." << std::endl;

  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Dirichlet boundary conditions
  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  if (!infiltrate_only_if_unfrozen_) {
    // Standard Neuman boundary conditions
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = bc->second;
    }
  } else {
    // Neumann boundary conditions that turn off if temp < freezing
    const Epetra_MultiVector& temp = *S_next_->GetFieldData("temperature")->ViewComponent("face");
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      if (temp[0][f] > 273.15) {
        bc_values_[f] = bc->second;
      } else {
        bc_values_[f] = 0.;
      }
    }
  }

  // seepage face -- pressure <= p_atm, outward mass flux >= 0
  const Epetra_MultiVector pressure = *S_next_->GetFieldData(key_)->ViewComponent("face");
  const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");
  for (bc=bc_seepage_->begin(); bc!=bc_seepage_->end(); ++bc) {
    int f = bc->first;
    if (pressure[0][f] < p_atm) {
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = bc->second;
    } else {
      bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      bc_values_[f] = p_atm;
    }
  }

  // surface coupling
  if (coupled_to_surface_via_head_) {
    // Face is Dirichlet with value of surface head
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& head = *S_next_->GetFieldData("surface_pressure")
        ->ViewComponent("cell",false);

    int ncells_surface = head.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to dirichlet
      bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      bc_values_[f] = head[0][c];
    }
  }

  // surface coupling
  if (coupled_to_surface_via_flux_) {
    // Face is Neumann with value of surface residual
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& flux = *S_next_->GetFieldData("surface_subsurface_flux")
        ->ViewComponent("cell",false);

    int ncells_surface = flux.MyLength();
    for (int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to Neumann
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = flux[0][c] / mesh_->face_area(f);
      // NOTE: flux[0][c] is in units of mols / s, where as Neumann BCs are in
      //       units of mols / s / A.  The right A must be chosen, as it is
      //       the subsurface mesh's face area, not the surface mesh's cell
      //       area.
    }
  }
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
Richards::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& pres) {
  Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
  int nfaces = pres_f.MyLength();
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET) {
      pres_f[0][f] = bc_values_[f];
    }
  }
};


bool Richards::modify_predictor(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "Modifying predictor:" << std::endl;

  bool changed(false);
  if (modify_predictor_bc_flux_ ||
      (modify_predictor_first_bc_flux_ && S_next_->cycle() == 0)) {
    changed |= ModifyPredictorFluxBCs_(h,u);
  }


  if (modify_predictor_wc_) {
    changed |= ModifyPredictorWC_(h,u);
  }

  if (modify_predictor_with_consistent_faces_) {
    changed |= ModifyPredictorConsistentFaces_(h,u);
  }
  return changed;
}

bool Richards::ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
            wrms_, &bc_markers_, &bc_values_));
  }

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  UpdatePermeabilityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData("numerical_rel_perm");
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  flux_predictor_->modify_predictor(h, u);
  changed_solution(); // mark the solution as changed, as modifying with
                      // consistent faces will then get the updated boundary
                      // conditions
  return true;
}

bool Richards::ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = getOSTab();
  std::cout << std::setprecision(15);
  std::cout << " old pressure top face = " << (*u->data())("face",500) << std::endl;

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "  modifications for consistent face pressures." << std::endl;
  CalculateConsistentFaces(u->data().ptr());
  std::cout << " new pressure top face = " << (*u->data())("face",500) << std::endl;
  return true;
}

bool Richards::ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u) {
  ASSERT(0);
  return false;

  // Teuchos::OSTab tab = getOSTab();
  // if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
  //   *out_ << "  modifications for consistent face pressures." << std::endl;

  // // assumes mostly incompressible liquid, and uses just ds/dt
  // Epetra_MultiVector& u_c = u->data()->ViewComponent("cell",false);

  // // pres at previous step
  // const Epetra_MultiVector& p1 = *S_inter_->GetFieldData("saturation_liquid")
  //   ->ViewComponent("cell",false);

  // // project saturation
  // double dt_next = S_next_->time() - S_inter_->time();
  // double dt_prev = S_inter_->time() - time_prev2_;

  // // -- get sat data
  // const Epetra_MultiVector& sat0 = *sat_prev2_;
  // const Epetra_MultiVector& sat1 = *S_inter_->GetFieldData("saturation_liquid")
  //     ->ViewComponent("cell",false);
  // Epetra_MultiVector& sat2 = *S_next_->GetFieldData("saturation_liquid",
  //         "saturation_liquid")->ViewComponent("cell",false);

  // // -- project
  // sat2 = sat0;
  // double dt_ratio = (dt_next + dt_prev) / dt_prev;
  // sat2.Update(dt_ratio, sat1, 1. - dt_ratio);

  // // determine pressure at new step
  // Epetra_MultiVector pres_wc(u_c);
  // int ncells = pres_wc.MyLength();
  // for (int c=0; c!=ncells; ++c) {
  //   // Query the WRM somehow?
  // }
}


void Richards::CalculateConsistentFacesForInfiltration_(
    const Teuchos::Ptr<CompositeVector>& u) {
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
            wrms_, &bc_markers_, &bc_values_));
  }

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  bool update = UpdatePermeabilityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  flux_predictor_->modify_predictor(u);
}

void Richards::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  // update the rel perm according to the scheme of choice
  changed_solution();
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  Teuchos::RCP<CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm", name_);

  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");

  std::cout << " Krel face 507 = " << (*rel_perm)("face",507) << std::endl;
  std::cout << " Pres cell 101 = " << (*u)("cell",101) << std::endl;


  // Update the preconditioner with darcy and gravity fluxes
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());

  // skip accumulation terms, they're not needed

  // Assemble
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // derive the consistent faces, involves a solve
  matrix_->UpdateConsistentFaceConstraints(u.ptr());
}

} // namespace
} // namespace
