/* -*-  mode++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Solves:

de/dt + q dot grad h = div Ke grad T + S?
------------------------------------------------------------------------- */

#include "advection.hh"
#include "FieldEvaluator.hh"
#include "energy_base.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Accumulation of energy term de/dt
// -------------------------------------------------------------
void EnergyBase::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the energy at both the old and new times.
  S_next_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_inter_.ptr(), name_);

  // get the energy at each time
  Teuchos::RCP<const CompositeVector> e1 = S_next_->GetFieldData(energy_key_);
  Teuchos::RCP<const CompositeVector> e0 = S_inter_->GetFieldData(energy_key_);

  // Update the residual with the accumulation of energy over the
  // timestep, on cells.
  g->ViewComponent("cell", false)
    ->Update(1.0/dt, *e1->ViewComponent("cell", false),
          -1.0/dt, *e0->ViewComponent("cell", false), 1.0);
};


// -------------------------------------------------------------
// Advective term for transport of enthalpy, q dot grad h.
// -------------------------------------------------------------
void EnergyBase::AddAdvection_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g, bool negate) {

  Teuchos::RCP<CompositeVector> field = advection_->field();
  field->PutScalar(0);

  // set the flux field
  // NOTE: fluxes are a MOLAR flux by choice of the flow pk, i.e.
  // [flux] =  mol/s

  // NOTE: this will be the eventual way to ensure it is up to date,
  // but there is no FieldEvaluator for darcy flux yet.  When there
  // is, we can take the evaluation out of Flow::commit_state(),
  // but for now we'll leave it there and assume it has been updated. --etc
  //  S->GetFieldEvaluator(flux_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(flux_key_);
  advection_->set_flux(flux);

  // put the advected quantity in cells
  S->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> enthalpy = S->GetFieldData(enthalpy_key_);
  *field->ViewComponent("cell", false) = *enthalpy->ViewComponent("cell", false);

  // put the boundary fluxes in faces for Dirichlet BCs.
  ApplyDirichletBCsToEnthalpy_(S.ptr(), field.ptr());

  // apply the advection operator and add to residual
  advection_->Apply(bc_flux_);

  Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);
  Teuchos::RCP<const CompositeVector> field_const(field);
  const Epetra_MultiVector& field_c =
      *field_const->ViewComponent("cell", false);

  int c_owned = g_c.MyLength();
  if (negate) {
    for (int c=0; c!=c_owned; ++c) {
      g_c[0][c] -= field_c[0][c];
    }
  } else {
    for (int c=0; c!=c_owned; ++c) {
      g_c[0][c] += field_c[0][c];
    }
  }
};


// -------------------------------------------------------------
// Diffusion term, div K grad T
// -------------------------------------------------------------
void EnergyBase::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g) {
  // update the thermal conductivity
  S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> conductivity =
    S->GetFieldData(conductivity_key_);

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(conductivity.ptr());
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

  // update the flux if needed
  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData(energy_flux_key_, name_);
    matrix_->DeriveFlux(*temp, flux.ptr());
  }

  // finish assembly of the stiffness matrix
  matrix_->CreateMFDrhsVectors();
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  matrix_->ComputeNegativeResidual(*temp, g);
};

} //namespace Energy
} //namespace Amanzi
