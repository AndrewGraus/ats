/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>

#include "ReconstructionCell.hh"
#include "OperatorDefs.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************* 
 * Routine takes a parallel overlapping vector C and returns a parallel
 * overlapping vector F(C).
 ****************************************************************** */
void Transport_PK::Functional(const double t, const Epetra_Vector& component, Epetra_Vector& f_component)
{
  // transport routines need an RCP pointer
  Teuchos::RCP<const Epetra_Vector> component_rcp(&component, false);

  Teuchos::ParameterList plist = tp_list_->sublist("reconstruction");
  lifting_->Init(component_rcp, plist);
  lifting_->Compute();
  Teuchos::RCP<CompositeVector> gradient = lifting_->gradient();

  // extract boundary conditions for the current component
  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost);

  for (int m = 0; m < bcs.size(); m++) {
    std::vector<int>& tcc_index = bcs[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        std::vector<int>& faces = bcs[m]->faces();
        std::vector<std::vector<double> >& values = bcs[m]->values();
        int nbfaces = faces.size();

        for (int n = 0; n < nbfaces; ++n) {
          int f = faces[n];
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[n][i];
        }
      }
    }
  }

  lifting_->InitLimiter(darcy_flux);
  lifting_->ApplyLimiter(bc_model, bc_value);

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  int f, c1, c2;
  double u, u1, u2, umin, umax, upwind_tcc, tcc_flux;

  f_component.PutScalar(0.0);
  for (int f = 0; f < nfaces_wghost; f++) {  // loop over master and slave faces
    c1 = (*upwind_cell_)[f];
    c2 = (*downwind_cell_)[f];

    if (c1 >= 0 && c2 >= 0) {
      u1 = component[c1];
      u2 = component[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = component[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = component[c2];
    }

    u = fabs((*darcy_flux)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    if (c1 >= 0 && c1 < ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = lifting_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;
    } else if (c1 >= 0 && c1 < ncells_owned && (c2 >= ncells_owned || c2 < 0)) {
      upwind_tcc = lifting_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
    } else if (c1 >= ncells_owned && c2 >= 0 && c2 < ncells_owned) {
      upwind_tcc = lifting_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  }

  // process external sources
  if (srcs.size() != 0) {
    ComputeAddSourceTerms(t, 1.0, srcs, f_component, current_component_, current_component_);
  }

  for (int c = 0; c < ncells_owned; c++) {  // calculate conservative quantatity
    double vol_phi_ws = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_start)[0][c];
    f_component[c] /= vol_phi_ws;
  }

  // BOUNDARY CONDITIONS for ADVECTION
  for (int m = 0; m < bcs.size(); m++) {
    std::vector<int>& tcc_index = bcs[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        std::vector<int>& faces = bcs[m]->faces();
        std::vector<std::vector<double> >& values = bcs[m]->values();
        int nbfaces = faces.size();

        for (int n = 0; n < nbfaces; ++n) {
          f = faces[n];
          c2 = (*downwind_cell_)[f];

          if (c2 >= 0 && f < nfaces_owned) {
            u = fabs((*darcy_flux)[0][f]);
            double vol_phi_ws = mesh_->cell_volume(c2) * (*phi)[0][c2] * (*ws_start)[0][c2];
            tcc_flux = u * values[n][i];
            f_component[c2] += tcc_flux / vol_phi_ws;
          }
        }
      }
    }
  }
}

}  // namespace Transport
}  // namespace Amanzi


