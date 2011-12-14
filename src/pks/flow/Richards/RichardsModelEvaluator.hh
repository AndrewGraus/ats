#ifndef RICHARDS_MODEL_EVALUATOR_HPP
#define RICHARDS_MODEL_EVALUATOR_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"

#include "BDF2_fnBase.hpp"
#include "RichardsProblem.hh"

namespace Amanzi {

class RichardsModelEvaluator : public BDF2::fnBase,
			       public Teuchos::VerboseObject<RichardsModelEvaluator> {
public:

  // Constructor
  RichardsModelEvaluator(Teuchos::RCP<RichardsProblem>& problem,
                         Teuchos::ParameterList& plist);
  // Initialization
  void initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params);

  void fun(const double t, const Epetra_Vector& u, const Epetra_Vector& udot, Epetra_Vector& f);
  void precon(const Epetra_Vector& u, Epetra_Vector& Pu);
  double enorm(const Epetra_Vector& u, const Epetra_Vector& du);
  void update_precon(const double t, const Epetra_Vector& up, const double h, int& errc);

  bool is_admissible(const Epetra_Vector& up);

private:
  Teuchos::RCP<RichardsProblem> problem_;

  // ML preconditioner for the Schur complement system.
  ML_Epetra::MultiLevelPreconditioner *MLprec;

  Teuchos::ParameterList& plist_;
  double atol, rtol;

};

} // close namespace Amanzi

#endif // RICHARDS_MODEL_EVALUATOR_HPP 
