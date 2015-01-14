/*
Author: Ethan Coon

Painter's permafrost model with freezing point depression.

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_fpd_smoothed_permafrost_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMFPDSmoothedPermafrostModel::saturations(double pc_liq, double pc_ice,
        double T, double (&sats)[3]) {
  if (pc_liq <= 0.) { // saturated
    sats[0] = 0.; // gas
    sats[1] = wrm_->saturation(pc_ice); // liquid
    sats[2] = 1.0 - sats[1];  // ice
    ASSERT(sats[2] >= 0.);
  } else if (pc_ice <= pc_liq) {
    sats[2] = 0.; // ice
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[0] = 1.0 - sats[1];  // gas
    ASSERT(sats[0] >= 0.);
  } else {
    sats[1] = wrm_->saturation(pc_ice);
    sats[2] = 1. - sats[1] / wrm_->saturation(pc_liq);
    sats[0] = 1. - sats[1] - sats[2];
    ASSERT(sats[2] >= 0.);
    ASSERT(sats[0] >= 0.);
  }
}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.; // gas
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas
  } else {
    dsats[1] = 0.;
    dsats[2] = wrm_->saturation(pc_ice) / std::pow(wrm_->saturation(pc_liq),2)
        * wrm_->d_saturation(pc_liq);
    dsats[0] = - dsats[2];
  }
}


void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.; // gas
    dsats[1] = wrm_->d_saturation(pc_ice); // liquid
    dsats[2] = - dsats[1];  // ice
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
  } else {
    dsats[1] = wrm_->d_saturation(pc_ice);
    dsats[2] = - dsats[1] / wrm_->saturation(pc_liq);
    dsats[0] = - dsats[1] - dsats[2];
  }

}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dtemperature(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  dsats[0] = 0.;
  dsats[1] = 0.;
  dsats[2] = 0.;
}


} // namespace FlowRelations
} // namespace Flow
} // namespace Amanzi
