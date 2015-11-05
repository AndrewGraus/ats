#include "passive_tracer.hh"

namespace Amanzi {
namespace Transport {

// Methods for the BDF integrator
// -- residual
void PassiveTracer::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
        Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);

  // pointer-copy temperature into states and update any auxilary data
  solution_to_state(u_old, S_inter_);
  solution_to_state(u_new, S_next_);

  // get access to the solution
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // update bcs
  for (int i=0; i != bcs_->size(); ++i) (*bcs_)[i]->Compute(t_new);

  // accumulation term
  AddAccumulation_(res);

  // advection term
  AddAdvection_(res);
};

// -- preconditioning (currently none)
int PassiveTracer::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
  return 0;
};

// computes a norm on u-du and returns the result
double PassiveTracer::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> conc_vec = u->Data()->ViewComponent(false);
  Teuchos::RCP<const Epetra_MultiVector> conc_dot_vec = du->Data()->ViewComponent(false);

  for (int lcv=0; lcv != conc_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*conc_dot_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*conc_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
  }

  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return enorm_val;
};

} // namespace
} // namespace
