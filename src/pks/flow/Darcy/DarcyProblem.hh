#ifndef __DARCYPROBLEM_H__
#define __DARCYPROBLEM_H__

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"
#include "MimeticHexLocal.hpp"
#include "MimeticHex.hpp"

namespace Amanzi
{

class DarcyMatvec; // forward declaration
class BoundaryFunction; // forward declaration

class DarcyProblem
{
public:

  DarcyProblem(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
	       Teuchos::ParameterList& darcy_plist);
  ~DarcyProblem();

  // initialization
  void InitializeProblem(Teuchos::ParameterList& list);

  // Set the constant value of fluid density.
  void SetFluidDensity(double rho);

  // Set the constant value of fluid viscosity.
  void SetFluidViscosity(double mu);

  // Set the gravitational acceleration (3-vector).
  void SetGravity(const double g[3]);

  // Sets a constant (scalar) permeability.
  void SetPermeability(double k);

  // Sets a spatially variable (scalar) permeability, one value per cell.
  //void SetPermeability(const std::vector<double> &k);
  void SetPermeability(const Epetra_Vector &k);

  // Assemble the problem
  void Assemble();

  void ComputeF (const Epetra_Vector &X, Epetra_Vector &F);

  const Epetra_Vector& RHS() const { return *rhs_; }

  const Epetra_Map& Map() const { return *dof_map_; }

  const Epetra_Map& CellMap(bool ghost=false) const { return mesh_->cell_map(ghost); }
  const Epetra_Map& FaceMap(bool ghost=false) const { return mesh_->face_map(ghost); }

  Epetra_Vector* CreateCellView(const Epetra_Vector&) const;
  Epetra_Vector* CreateFaceView(const Epetra_Vector&) const;

  const Epetra_Comm& Comm() const { return *(mesh_->get_comm()); }

  Epetra_Operator& Precon() const { return *precon_; }

  Epetra_Operator& Matvec() const { return *matvec_; }

  DiffusionMatrix& Matrix() const { return *D_; }

  void DeriveDarcyFlux(const Epetra_Vector &P, Epetra_Vector &F, double &l1_error) const;
  void DeriveDarcyVelocity(const Epetra_Vector &X, Epetra_MultiVector &Q) const;

  void GetFluidDensity(double &rho) const { rho = rho_; }
  void GetFluidViscosity(double &mu) const { mu = mu_; }
  void GetGravity(double g[]) const { for(int i = 0; i < 3; ++i) g[i] = g_[i]; }

private:

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Epetra_Map> dof_map_;
  Teuchos::RCP<Epetra_Import> face_importer_;
  Teuchos::RCP<Epetra_Import> cell_importer_;

  double rho_;  // constant fluid density
  double mu_;   // constant fluid viscosity
  double g_[3]; // gravitational acceleration
  double gravity_;  // gravitational acceleration (positive coef, directed in -z direction)
  std::vector<double> k_; // spatially variable permeability

  Teuchos::RCP<BoundaryFunction> bc_press_;  // Pressure Dirichlet conditions, excluding static head
  Teuchos::RCP<BoundaryFunction> bc_head_;   // Static pressure head conditions; also Dirichlet-type
  Teuchos::RCP<BoundaryFunction> bc_flux_;   // Outward mass flux conditions

  std::vector<MimeticHexLocal>  MD_;
  MimeticHex *md_; // evolving replacement for MD

  Epetra_Vector *rhs_;

  DiffusionPrecon *precon_;
  Epetra_Operator *matvec_;

  Teuchos::RCP<DiffusionMatrix> D_;

private:  // Auxillary functions

  Teuchos::RCP<Epetra_Map> create_dof_map_(const Epetra_Map&, const Epetra_Map&) const;
  void validate_boundary_conditions() const;
  DiffusionMatrix* create_diff_matrix_(Teuchos::RCP<AmanziMesh::Mesh>&) const;
  void init_mimetic_disc_(Teuchos::RCP<AmanziMesh::Mesh>&, std::vector<MimeticHexLocal>&) const;
};

} // close namespace Amanzi

#endif
