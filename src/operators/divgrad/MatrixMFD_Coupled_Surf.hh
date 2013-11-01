/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  Class to form the 2x n_lambda system for coupling energy and flow on both
  the surface and subsurface.

*/

#ifndef OPERATORS_MATRIX_COUPLED_MFD_SURF_HH_
#define OPERATORS_MATRIX_COUPLED_MFD_SURF_HH_


#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD_Coupled.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Coupled_Surf : public MatrixMFD_Coupled {

 public:
  MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh> mesh);
  MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other);

  virtual void ComputeSchurComplement();

  void SetOffDiagonals(const Teuchos::RCP<const Epetra_MultiVector>& Ccc,
                       const Teuchos::RCP<const Epetra_MultiVector>& Dcc,
                       const Teuchos::RCP<const Epetra_MultiVector>& Ccc_surf=Teuchos::null,
                       const Teuchos::RCP<const Epetra_MultiVector>& Dcc_surf=Teuchos::null,
                       double scaling=1.) {
    scaling_ = scaling;
    MatrixMFD_Coupled::SetOffDiagonals(Ccc,Dcc,scaling);
    Ccc_surf_ = Ccc_surf;
    Dcc_surf_ = Dcc_surf;
  }


  void SetSurfaceOperators(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           const Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
    surface_A_ = surface_A;
    surface_B_ = surface_B;
    ASSERT(surface_A_->Mesh() == surface_B_->Mesh());
    surface_mesh_ = surface_A_->Mesh();
  }

  void GetSurfaceOperators(Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
    surface_A = surface_A_;
    surface_B = surface_B_;
  }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_B_;

  Teuchos::RCP<const Epetra_MultiVector> Ccc_surf_;
  Teuchos::RCP<const Epetra_MultiVector> Dcc_surf_;

  double scaling_;
};


} // namespace
} // namespace


#endif
