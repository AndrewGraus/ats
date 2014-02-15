#include "mpc_surface_subsurface_helpers.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf,
                        const Teuchos::Ptr<CompositeVector>& sub) {
  const Epetra_MultiVector& surf_c = *surf.ViewComponent("cell",false);
  Epetra_MultiVector& sub_f = *sub->ViewComponent("face",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    sub_f[0][f] = surf_c[0][sc];
  }
}

void
CopySubsurfaceToSurface(const CompositeVector& sub,
                        const Teuchos::Ptr<CompositeVector>& surf) {
  const Epetra_MultiVector& sub_f = *sub.ViewComponent("face",false);
  Epetra_MultiVector& surf_c = *surf->ViewComponent("cell",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    surf_c[0][sc] = sub_f[0][f];
  }
}

void
MergeSubsurfaceAndSurfacePressure(const Teuchos::Ptr<CompositeVector>& sub_p,
        const Teuchos::Ptr<CompositeVector>& surf_p,
        const CompositeVector& kr_surf) {
  Epetra_MultiVector& sub_p_f = *sub_p->ViewComponent("face",false);
  Epetra_MultiVector& surf_p_c = *surf_p->ViewComponent("cell",false);
  const Epetra_MultiVector& kr_c = *kr_surf.ViewComponent("cell",false);
  double p_atm = 101325.;

  for (unsigned int sc=0; sc!=surf_p_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf_p->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    //    if (kr_c[0][sc] > 0.) {
    if (surf_p_c[0][sc] < 101325.) {
      surf_p_c[0][sc] = sub_p_f[0][f];
    } else {
      sub_p_f[0][f] = surf_p_c[0][sc];
    }
  }
}


} // namespace
