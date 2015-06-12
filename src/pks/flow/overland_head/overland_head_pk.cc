/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_MultiVector.h"

#include "flow_bc_factory.hh"
#include "Mesh.hh"
#include "Point.hh"
#include "Op.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"
#include "independent_variable_field_evaluator.hh"

#include "upwind_potential_difference.hh"
#include "upwind_total_flux.hh"
#include "upwind_cell_centered.hh"
#include "pres_elev_evaluator.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "overland_head_water_content_evaluator.hh"
#include "height_model.hh"
#include "height_evaluator.hh"
#include "overland_source_from_subsurface_flux_evaluator.hh"

#include "OperatorDiffusionFactory.hh"

#include "overland_head.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0

OverlandHeadFlow::OverlandHeadFlow(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist, FElist, solution),
    PKPhysicalBDFBase(plist, FElist, solution),
    standalone_mode_(false),
    is_source_term_(false),
    coupled_to_subsurface_via_head_(false),
    coupled_to_subsurface_via_flux_(false),
    perm_update_required_(true),
    update_flux_(UPDATE_FLUX_ITERATION),
    full_jacobian_(false),
    niter_(0),
    source_only_if_unfrozen_(false),
    precon_used_(true)
{
  plist_->set("primary variable key", "surface_pressure");
  plist_->set("conserved quantity key", "surface_water_content");
  plist_->set("domain name", "surface");
  
  // clone the ponded_depth parameter list for ponded_depth bar
  Teuchos::ParameterList& pd_list = FElist.sublist("ponded_depth");
  Teuchos::ParameterList pdbar_list(pd_list);
  pdbar_list.set("ponded depth bar", true);
  pdbar_list.set("height key", "ponded_depth_bar");
  FElist.set("ponded_depth_bar", pdbar_list);

  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .01 * 55000.); // h * nl

}


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
void OverlandHeadFlow::setup(const Teuchos::Ptr<State>& S) {
  // set up the meshes
  if (!S->HasMesh("surface")) {
    Teuchos::RCP<const AmanziMesh::Mesh> domain = S->GetMesh();
    //    ASSERT(domain->space_dimension() == 2);
    standalone_mode_ = true;
    S->AliasMesh("domain", "surface");
  } else {
    standalone_mode_ = false;
  }

  // -- water content
  S->RequireField("surface_water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList& wc_plist =
      plist_->sublist("overland water content evaluator");
  Teuchos::RCP<FlowRelations::OverlandHeadWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadWaterContentEvaluator(wc_plist));
  S->SetFieldEvaluator("surface_water_content", wc_evaluator);

  PKPhysicalBDFBase::setup(S);
  SetupOverlandFlow_(S);
  SetupPhysicalEvaluators_(S);
}


void OverlandHeadFlow::SetupOverlandFlow_(const Teuchos::Ptr<State>& S) {
  // Require fields and evaluators for those fields.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2, 1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  
  
  std::vector<AmanziMesh::Entity_kind> locations1(1);
  std::vector<std::string> names1(1);
  std::vector<int> num_dofs1(1,1);
  locations1[0] = AmanziMesh::CELL;
  names1[0] = "cell";
  
  if (!(plist_->sublist("Diffusion").get<bool>("TPFA use cells only", false)) ){
    // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
    S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
    ->SetComponents(names2, locations2, num_dofs2);
  }
  else {
    // -- primary variable: pressure on both cells and boundary faces, ghosted, with 1 dof
    S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
    ->SetComponents(names1, locations1, num_dofs1);
  }
  
// #if DEBUG_RES_FLAG
//   // -- residuals of various iterations for debugging
//   for (int i=1; i!=23; ++i) {
//     std::stringstream namestream;
//     namestream << "flow_residual_" << i;
//     std::stringstream solnstream;
//     solnstream << "flow_solution_" << i;
//     S->RequireField(namestream.str(), name_)->SetMesh(mesh_)->SetGhosted()
//                     ->SetComponents(names2, locations2, num_dofs2);
//     S->RequireField(solnstream.str(), name_)->SetMesh(mesh_)->SetGhosted()
//                     ->SetComponents(names2, locations2, num_dofs2);
//   }
// #endif

  // -- cell volume and evaluator
  S->RequireFieldEvaluator("surface_cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // -- coupling to subsurface
  coupled_to_subsurface_via_flux_ =
      plist_->get<bool>("coupled to subsurface via flux", false);
  coupled_to_subsurface_via_head_ =
      plist_->get<bool>("coupled to subsurface via head", false);
  ASSERT(!(coupled_to_subsurface_via_flux_ && coupled_to_subsurface_via_head_));

  if (coupled_to_subsurface_via_head_) {
    // -- source term from subsurface, filled in by evaluator,
    //    which picks the fluxes from "darcy_flux" field.
    S->RequireFieldEvaluator("surface_subsurface_flux");
    S->RequireField("surface_subsurface_flux")
        ->SetMesh(mesh_)->SetComponents(names1, locations1, num_dofs1);
  }

  // Create the upwinding method.
  S->RequireField("upwind_overland_conductivity", name_)->SetMesh(mesh_)
      ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->GetField("upwind_overland_conductivity",name_)->set_io_vis(false);

  S->RequireField("dupwind_overland_conductivity_dponded_depth", name_)->SetMesh(mesh_)
      ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  S->GetField("dupwind_overland_conductivity_dponded_depth",name_)->set_io_vis(false);

  upwind_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
          "overland_conductivity", "upwind_overland_conductivity",
          "surface_flux_direction", 1.e-8));
    
  upwinding_dkdp_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
          "doverland_conductivity_dponded_depth",
          "dupwind_overland_conductivity_dponded_depth",
          "surface_flux_direction", 1.e-8));

  // -- owned secondary variables, no evaluator used
  S->RequireField("surface_flux_direction", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("surface_flux", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_head_ = bc_factory.CreateHead();
  bc_zero_gradient_ = bc_factory.CreateZeroGradient();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_seepage_head_ = bc_factory.CreateWithFunction("seepage face head", "boundary head");
  bc_seepage_pressure_ = bc_factory.CreateWithFunction("seepage face pressure", "boundary pressure");
  ASSERT(!bc_plist.isParameter("seepage face")); // old style!

  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);
  std::vector<double> mixed;
  bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_markers_, bc_values_, mixed));

  // operator for the diffusion terms: must use ScaledConstraint version
  Teuchos::ParameterList& mfd_plist = plist_->sublist("Diffusion");
  mfd_plist.set("scaled constraint equation", true);
  Operators::OperatorDiffusionFactory opfactory;
  matrix_diff_ = opfactory.Create(mesh_, bc_, mfd_plist);
  matrix_diff_->Setup(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // operator for flux directions -- this should be removed eventually
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_diff_list.set("scaled constraint equation", true);
  face_matrix_diff_ = opfactory.Create(mesh_, bc_, face_diff_list);
  face_matrix_diff_->Setup(Teuchos::null);
  face_matrix_diff_->Setup(Teuchos::null, Teuchos::null);
  face_matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  
  // diffusion operator for the preconditioner
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("Diffusion PC");
  mfd_pc_plist.set("scaled constraint equation", true);
  preconditioner_diff_ = opfactory.Create(mesh_, bc_, mfd_pc_plist);
  preconditioner_diff_->Setup(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();

  // accumulation operator for the preconditioenr
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("Accumulation PC");
  acc_pc_plist.set("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(acc_pc_plist, preconditioner_));

  precon_used_ = plist_->isSublist("preconditioner");
  if (precon_used_) {
    preconditioner_->SymbolicAssembleMatrix();
  }

  // how often to update the fluxes?
  update_flux_ = UPDATE_FLUX_ITERATION;
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void OverlandHeadFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2, 1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  
  
  std::vector<AmanziMesh::Entity_kind> locations1(1);
  std::vector<std::string> names1(1);
  std::vector<int> num_dofs1(1,1);
  locations1[0] = AmanziMesh::CELL;
  names1[0] = "cell";
  
  // -- evaluator for surface geometry.
  S->RequireField("elevation")->SetMesh(S->GetMesh("surface"))->SetGhosted()
  ->AddComponents(names2, locations2, num_dofs2);

  S->RequireField("slope_magnitude")->SetMesh(S->GetMesh("surface"))
      ->AddComponents(names1, locations1, num_dofs1);

  Teuchos::RCP<FlowRelations::ElevationEvaluator> elev_evaluator;
  if (standalone_mode_) {
    ASSERT(plist_->isSublist("elevation evaluator"));
    Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::StandaloneElevationEvaluator(elev_plist));
  } else {
    Teuchos::ParameterList elev_plist = plist_->sublist("elevation evaluator");
    elev_evaluator = Teuchos::rcp(new FlowRelations::MeshedElevationEvaluator(elev_plist));
  }
  S->SetFieldEvaluator("elevation", elev_evaluator);
  S->SetFieldEvaluator("slope_magnitude", elev_evaluator);

  // -- evaluator for potential field, h + z
  if (!(plist_->sublist("Diffusion").get<bool>("TPFA use cells only", false)) ){
    S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
    ->SetComponents(names2, locations2, num_dofs2);
  }
  else {
    S->RequireField("pres_elev")->SetMesh(S->GetMesh("surface"))->SetGhosted()
    ->SetComponents(names1, locations1, num_dofs1);
  }
  Teuchos::ParameterList pres_elev_plist = plist_->sublist("potential evaluator");
  Teuchos::RCP<FlowRelations::PresElevEvaluator> pres_elev_eval =
      Teuchos::rcp(new FlowRelations::PresElevEvaluator(pres_elev_plist));
  S->SetFieldEvaluator("pres_elev", pres_elev_eval);

  // -- evaluator for source term
  is_source_term_ = plist_->get<bool>("source term");
  if (is_source_term_) {
    source_in_meters_ = plist_->get<bool>("mass source in meters", true);
    
    // source term itself [m/s]
    mass_source_key_ = plist_->get<std::string>("mass source key", "surface_mass_source");
    S->RequireField(mass_source_key_)->SetMesh(mesh_)
        ->AddComponents(names1, locations1, num_dofs1);
    S->RequireFieldEvaluator(mass_source_key_);

    // density of incoming water [mol/m^3]
    S->RequireField("surface_source_molar_density")->SetMesh(mesh_)
        ->AddComponents(names1, locations1, num_dofs1);
    S->RequireFieldEvaluator("surface_source_molar_density");
  }

  // -- water content bar (can be negative)
  S->RequireField("surface_water_content_bar")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names1, locations1, num_dofs1);
  Teuchos::ParameterList& wc_plist =
      plist_->sublist("overland water content evaluator");
  Teuchos::ParameterList wcbar_plist(wc_plist);
  wcbar_plist.set<bool>("water content bar", true);
  Teuchos::RCP<FlowRelations::OverlandHeadWaterContentEvaluator> wc_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandHeadWaterContentEvaluator(wcbar_plist));
  S->SetFieldEvaluator("surface_water_content_bar", wc_evaluator);

  // -- ponded depth
  if (!(plist_->sublist("Diffusion").get<bool>("TPFA use cells only", false)) ){
    S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
    ->SetComponents(names2, locations2, num_dofs2);
  }
  else {
    S->RequireField("ponded_depth")->SetMesh(mesh_)->SetGhosted()
    ->SetComponents(names1, locations1, num_dofs1);
  }
  S->RequireFieldEvaluator("ponded_depth");

  // -- ponded depth bar
  S->RequireField("ponded_depth_bar")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names1, locations1, num_dofs1);
  S->RequireFieldEvaluator("ponded_depth_bar");

  // -- effective accumulation ponded depth (smoothing of derivatives as h --> 0)
  smoothed_ponded_accumulation_ = plist_->get<bool>("smooth ponded accumulation",false);
  if (smoothed_ponded_accumulation_) {
    S->RequireField("smoothed_ponded_depth")->SetMesh(mesh_)
      ->AddComponents(names1, locations1, num_dofs1);
    S->RequireFieldEvaluator("smoothed_ponded_depth");
  }

  // -- conductivity evaluator
  S->RequireField("overland_conductivity")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names1, locations1, num_dofs1);
  ASSERT(plist_->isSublist("overland conductivity evaluator"));
  Teuchos::ParameterList cond_plist = plist_->sublist("overland conductivity evaluator");
  Teuchos::RCP<FlowRelations::OverlandConductivityEvaluator> cond_evaluator =
      Teuchos::rcp(new FlowRelations::OverlandConductivityEvaluator(cond_plist));
  S->SetFieldEvaluator("overland_conductivity", cond_evaluator);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void OverlandHeadFlow::initialize(const Teuchos::Ptr<State>& S) {
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

  Teuchos::RCP<CompositeVector> head_cv = S->GetFieldData(key_, name_);

  // initial condition is tricky
  if (!S->GetField(key_)->initialized()) {
    if (!plist_->isSublist("initial condition")) {
      std::stringstream messagestream;
      messagestream << name_ << " has no initial condition parameter list.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
    head_cv->PutScalar(0.);
  }

  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  if (!S->GetField(key_)->initialized()) {
    // -- set the cell initial condition if it is taken from the subsurface
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
      Epetra_MultiVector& head = *head_cv->ViewComponent("cell",false);
      const Epetra_MultiVector& pres = *S->GetFieldData("pressure")
        ->ViewComponent("face",false);

      unsigned int ncells_surface = mesh_->num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
      for (unsigned int c=0; c!=ncells_surface; ++c) {
        // -- get the surface cell's equivalent subsurface face and neighboring cell
        AmanziMesh::Entity_ID f =
          mesh_->entity_get_parent(AmanziMesh::CELL, c);
        head[0][c] = pres[0][f];
      }

      // mark as initialized
      if (ic_plist.get<bool>("initialize surface head from subsurface",false)) {
        S->GetField(key_,name_)->set_initialized();
      }
    }
  }

  // Initialize BC data structures
  unsigned int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);
  std::vector<double> mixed;
  bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_markers_, bc_values_, mixed));

  // Initialize BC values
  bc_head_->Compute(S->time());
  bc_zero_gradient_->Compute(S->time());
  bc_flux_->Compute(S->time());
  bc_seepage_head_->Compute(S->time());
  bc_seepage_pressure_->Compute(S->time());

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("upwind_overland_conductivity",name_)->PutScalar(1.0);
  S->GetField("upwind_overland_conductivity",name_)->set_initialized();
  S->GetFieldData("dupwind_overland_conductivity_dponded_depth",name_)->PutScalar(1.0);
  S->GetField("dupwind_overland_conductivity_dponded_depth",name_)->set_initialized();
  S->GetField("surface_flux", name_)->set_initialized();
  S->GetFieldData("surface_flux_direction", name_)->PutScalar(0.);
  S->GetField("surface_flux_direction", name_)->set_initialized();
  //  S->GetField("surface_velocity", name_)->set_initialized();
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::commit_state(double dt, const Teuchos::RCP<State>& S) {
  niter_ = 0;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;

  PKPhysicalBDFBase::commit_state(dt, S);

  // update boundary conditions
  bc_head_->Compute(S->time());
  bc_flux_->Compute(S->time());
  bc_seepage_head_->Compute(S->time());
  bc_seepage_pressure_->Compute(S->time());
  UpdateBoundaryConditions_(S.ptr());

  // Update flux if rel perm or h + Z has changed.
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S.ptr(), name_);

  // update the stiffness matrix with the new rel perm
  Teuchos::RCP<const CompositeVector> conductivity =
      S->GetFieldData("upwind_overland_conductivity");
  matrix_->Init();
  matrix_diff_->Setup(conductivity, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // Patch up BCs for zero-gradient
  FixBCsForOperator_(S.ptr());
  
  // derive the fluxes
  Teuchos::RCP<const CompositeVector> potential = S->GetFieldData("pres_elev");
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData("surface_flux", name_);
  matrix_diff_->UpdateFlux(*potential, *flux);

};


// -----------------------------------------------------------------------------
// Update diagnostics -- used prior to vis.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::calculate_diagnostics(const Teuchos::RCP<State>& S) {};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool OverlandHeadFlow::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S, name_);
  update_perm |= S->GetFieldEvaluator("pres_elev")->HasFieldChanged(S, name_);
  update_perm |= perm_update_required_;

  if (update_perm) {
    // Update the perm only if needed.
    perm_update_required_ = false;

    // get upwind conductivity data
    Teuchos::RCP<CompositeVector> uw_cond =
        S->GetFieldData("upwind_overland_conductivity", name_);

    // update the direction of the flux -- note this is NOT the flux
    Teuchos::RCP<CompositeVector> flux_dir =
        S->GetFieldData("surface_flux_direction", name_);
    Teuchos::RCP<const CompositeVector> pres_elev = S->GetFieldData("pres_elev");
    face_matrix_diff_->UpdateFlux(*pres_elev, *flux_dir);

    // get conductivity data
    Teuchos::RCP<const CompositeVector> cond = S->GetFieldData("overland_conductivity");
    const Epetra_MultiVector& cond_c = *cond->ViewComponent("cell",false);

    // place internal cell's value on faces -- this should be fixed to be the boundary data
    { // place boundary_faces on faces
      Epetra_MultiVector& uw_cond_f = *uw_cond->ViewComponent("face",false);

      AmanziMesh::Entity_ID_List cells;
      int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        if (cells.size() == 1) {
          int c = cells[0];
          uw_cond_f[0][f] = cond_c[0][c];
        }
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_->Update(S);
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Derivatives of the overland conductivity, upwinded.
// -----------------------------------------------------------------------------
bool OverlandHeadFlow::UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  // bool update_perm = S->GetFieldEvaluator("overland_conductivity")
  //     ->HasFieldDerivativeChanged(S, name_, key_);
  bool update_perm = S->GetFieldEvaluator("overland_conductivity")
      ->HasFieldDerivativeChanged(S, name_, "ponded_depth");

  if (update_perm) {
    // get upwind conductivity data
    Teuchos::RCP<CompositeVector> duw_cond =
        S->GetFieldData("dupwind_overland_conductivity_dponded_depth", name_);

    // get conductivity data
    Teuchos::RCP<const CompositeVector> dcond = S->GetFieldData("doverland_conductivity_dponded_depth");
    const Epetra_MultiVector& dcond_c = *dcond->ViewComponent("cell",false);

    { // place boundary_faces on faces
      Epetra_MultiVector& duw_cond_f = *duw_cond->ViewComponent("face",false);

      AmanziMesh::Entity_ID_List cells;
      int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
      for (int f=0; f!=nfaces_owned; ++f) {
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        if (cells.size() == 1) {
          int c = cells[0];
          duw_cond_f[0][f] = dcond_c[0][c];
        }
      }
    }

    // Then upwind.  This overwrites the boundary if upwinding says so.
    upwinding_dkdp_->Update(S);
  }

  if (update_perm && vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << " TRUE." << std::endl;
  return update_perm;
}


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void OverlandHeadFlow::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
      ->ViewComponent("face",false);

  // initialize all as null
  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_[n] = 0.0;
  }

  // Head BCs are standard Dirichlet, plus elevation
  for (Functions::BoundaryFunction::Iterator bc=bc_head_->begin();
       bc!=bc_head_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values_[f] = bc->second + elevation[0][f];
  }

  // Standard Neumann data for flux
  for (Functions::BoundaryFunction::Iterator bc=bc_flux_->begin();
       bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_[f] = bc->second;
  }

  // zero gradient: grad h = 0 implies that q = -k grad z
  // -- cannot be done yet as rel perm update is done after this and is needed.
  // -- Instead zero gradient BCs are done in FixBCs methods.

  // Seepage face head boundary condition
  if (bc_seepage_head_->size() > 0) {
    S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_c = *S->GetFieldData("ponded_depth")->ViewComponent("cell");
    const Epetra_MultiVector& elevation_c = *S->GetFieldData("elevation")->ViewComponent("cell");

    for (Functions::BoundaryFunction::Iterator bc = bc_seepage_head_->begin(); 
         bc != bc_seepage_head_->end(); ++bc) {
      int f = bc->first;
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];

      double hz_f = bc->second + elevation[0][f];
      double hz_c = h_c[0][c] + elevation_c[0][f];

      if (hz_f >= hz_c) {
        bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_values_[f] = 0.0;
      } else {
        bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_values_[f] = hz_f;
      }
    }
  }

  // Seepage face pressure boundary condition
  if (bc_seepage_head_->size() > 0) {
    S->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S.ptr(), name_);

    const Epetra_MultiVector& h_cells = *S->GetFieldData("ponded_depth")->ViewComponent("cell");
    const Epetra_MultiVector& elevation_cells = *S->GetFieldData("elevation")->ViewComponent("cell");
    const Epetra_MultiVector& rho_l = *S->GetFieldData("surface_mass_density_liquid")->ViewComponent("cell");
    double gz = -(*S->GetConstantVectorData("gravity"))[2];
    const double& p_atm = *S->GetScalarData("atmospheric_pressure");

    if (S->HasFieldEvaluator("surface_mass_density_ice")) {
      // thermal model of height
      const Epetra_MultiVector& eta = *S->GetFieldData("unfrozen_fraction")->ViewComponent("cell");
      const Epetra_MultiVector& rho_i = *S->GetFieldData("surface_mass_density_ice")->ViewComponent("cell");

      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_pressure_->begin(); 
           bc != bc_seepage_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / ((eta[0][c]*rho_l[0][c] + (1.-eta[0][c])*rho_i[0][c]) * gz);
        double dz = elevation_cells[0][c] - elevation[0][f];

        if (h_cells[0][c] + dz < h0) {
          bc_markers_[f] = Operators::OPERATOR_BC_NONE;
          bc_values_[f] = 0.0;
        } else {
          bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_values_[f] = h0 + elevation[0][f];
        }
      }

    } else {
      // non-thermal model
      for (Functions::BoundaryFunction::Iterator bc = bc_seepage_pressure_->begin(); 
           bc != bc_seepage_pressure_->end(); ++bc) {
        int f = bc->first;
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        int c = cells[0];

        double p0 = bc->second > p_atm ? bc->second : p_atm;
        double h0 = (p0 - p_atm) / (rho_l[0][c] * gz);
        double dz = elevation_cells[0][c] - elevation[0][f];

        if (h_cells[0][c] + dz < h0) {
          bc_markers_[f] = Operators::OPERATOR_BC_NONE;
          bc_values_[f] = 0.0;
        } else {
          bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_values_[f] = h0 + elevation[0][f];
        }
      }
    }
  }



  // mark all remaining boundary conditions as zero flux conditions
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_values_[f] = 0.0;
      }
    }
  }
  
}


void OverlandHeadFlow::FixBCsForOperator_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the Operator." << std::endl;

  // Now we can safely calculate q = -k grad z for zero-gradient problems
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData("elevation");
  elev->ScatterMasterToGhosted();
  const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face",false);
  const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell",false);

  std::vector<WhetStone::DenseMatrix>& Aff =
      matrix_diff_->local_matrices()->matrices;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
       bc!=bc_zero_gradient_->end(); ++bc) {

    int f = bc->first;

    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c = cells[0];

    if (f < nfaces_owned) {
      double dp = elevation_f[0][f] - elevation_c[0][c];
      double bc_val = -dp * Aff[f](0,0);

      bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_values_[f] = bc_val / mesh_->face_area(f);
    }
  }
};


void OverlandHeadFlow::FixBCsForPrecon_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "    Tweaking BCs for the PC." << std::endl;

  // // Attempt of a hack to deal with zero rel perm
  // double eps = 1.e-30;
  // Teuchos::RCP<CompositeVector> relperm =
  //     S->GetFieldData("upwind_overland_conductivity", name_);
  // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
  //   if ((*relperm)("face",f) < eps) {
  //     if (bc_markers_[f] == Operators::OPERATOR_BC_NEUMANN) {
  //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     } else if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
  //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
  //     }
  //   }
  // }
};

// void OverlandHeadFlow::FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S) {
//   Teuchos::OSTab tab = vo_->getOSTab();
//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "    Tweaking BCs for calculation of consistent faces." << std::endl;

//   // // If the rel perm is 0, the face value drops out and is unconstrained.
//   // // Therefore we set it to Dirichlet to eliminate it from the system.
//   // double eps = 1.e-30;
//   // const Epetra_MultiVector& elevation = *S->GetFieldData("elevation")
//   //     ->ViewComponent("face",false);
//   // Teuchos::RCP<CompositeVector> relperm =
//   //     S->GetFieldData("upwind_overland_conductivity", name_);

//   // for (unsigned int f=0; f!=relperm->size("face"); ++f) {
//   //   if ((*relperm)("face",f) < eps) {
//   //     if (bc_markers_[f] == Operators::OPERATOR_BC_NEUMANN) {
//   //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
//   //       bc_values_[f] =  elevation[0][f];
//   //     } else if (bc_markers_[f] == Operators::OPERATOR_BC_NONE) {
//   //       bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
//   //       bc_values_[f] =  elevation[0][f];
//   //     }
//   //   }
//   // }

//   // Now we can safely calculate q = -k grad z for zero-gradient problems
//   Teuchos::RCP<const CompositeVector> elev = S->GetFieldData("elevation");
//   elev->ScatterMasterToGhosted();
//   const Epetra_MultiVector& elevation_f = *elev->ViewComponent("face",false);
//   const Epetra_MultiVector& elevation_c = *elev->ViewComponent("cell",false);

//   std::vector<WhetStone::DenseMatrix>& Aff_cells =
//       matrix_diff_->local_matrices()->matrices;
//   Epetra_MultiVector& rhs_f = *matrix_->rhs()->ViewComponent("face",false);

//   int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
//   for (Functions::BoundaryFunction::Iterator bc=bc_zero_gradient_->begin();
//        bc!=bc_zero_gradient_->end(); ++bc) {
//     int f = bc->first;

//     AmanziMesh::Entity_ID_List cells;
//     mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
//     ASSERT(cells.size() == 1);
//     AmanziMesh::Entity_ID c = cells[0];

//     if (c < ncells_owned) {
//       AmanziMesh::Entity_ID_List faces;
//       mesh_->cell_get_faces(c, &faces);

//       std::vector<double> dp(faces.size());
//       for (unsigned int n=0; n!=faces.size(); ++n) {
//         dp[n] = elevation_f[0][faces[n]] - elevation_c[0][c];
//       }
//       unsigned int my_n = std::find(faces.begin(), faces.end(), f) - faces.begin();
//       ASSERT(my_n !=faces.size());

//       double bc_val = 0.;
//       for (unsigned int m=0; m!=faces.size(); ++m) {
//         bc_val -= Aff_cells[c](my_n,m) * dp[m];
//       }

//       // Apply the BC to the matrix
//       rhs_f[0][f] -= bc_val;
//     }
//   }
// };


/* ******************************************************************
 * Add a boundary marker to owned faces.
 ****************************************************************** */
void OverlandHeadFlow::ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& pres) {
  Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
  unsigned int nfaces = pres_f.MyLength();
  for (unsigned int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
      pres_f[0][f] = bc_values_[f];
    }
  }
};


bool OverlandHeadFlow::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                                       Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  // if (modify_predictor_with_consistent_faces_) {
  //   if (vo_->os_OK(Teuchos::VERB_EXTREME))
  //     *vo_->os() << "  modifications for consistent face pressures." << std::endl;
  //   CalculateConsistentFaces(u->Data().ptr());
  //   return true;
  // }

  return false;
};


void OverlandHeadFlow::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // // VerboseObject stuff.
  // Teuchos::OSTab tab = vo_->getOSTab();

  // // update the rel perm according to the scheme of choice
  // ChangedSolution();
  // UpdatePermeabilityData_(S_next_.ptr());

  // // update boundary conditions
  // bc_head_->Compute(S_next_->time());
  // bc_flux_->Compute(S_next_->time());
  // bc_seepage_head_->Compute(S_next_->time());
  // bc_seepage_pressure_->Compute(S_next_->time());
  // UpdateBoundaryConditions_(S_next_.ptr());

  // // update the stiffness matrix
  // Teuchos::RCP<const CompositeVector> cond =
  //   S_next_->GetFieldData("upwind_overland_conductivity", name_);
  // matrix_->CreateMFDstiffnessMatrices(cond.ptr());
  // matrix_->CreateMFDrhsVectors();

  // // Patch up BCs in the case of zero conductivity
  // FixBCsForConsistentFaces_(S_next_.ptr());

  // // Grab needed data.
  // S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
  // Teuchos::RCP<CompositeVector> pres_elev = S_next_->GetFieldData("pres_elev","pres_elev");

  // // Update the preconditioner with darcy and gravity fluxes
  // // skip accumulation terms, they're not needed
  // // Assemble
  // matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  // // derive the consistent faces, involves a solve
  // matrix_->UpdateConsistentFaceConstraints(pres_elev.ptr());

  // // back out heights from pres_elev
  // const Epetra_MultiVector& elevation = *S_next_->GetFieldData("elevation")
  //     ->ViewComponent("face",false);
  // u->ViewComponent("face",false)->Update(1., *pres_elev->ViewComponent("face",false),
  //         -1., elevation, 0.);
}

} // namespace
} // namespace

