/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_OVERLAND_CONDUCTIVITY_MODEL_

namespace Amanzi {
namespace Flow {

class OverlandConductivityModel {
 public:
  virtual ~OverlandConductivityModel() = default;
  
  virtual double Conductivity(double depth, double slope, double coef) = 0;
  virtual double DConductivityDDepth(double depth, double slope, double coef) = 0;

  //Add for the subgrid model -- Not pure virtual
  virtual double Conductivity(double depth, double slope, double coef, double d, double frac, double beta) {
    return Conductivity(depth,slope,coef);
  }
  virtual double DConductivityDDepth(double depth, double slope, double coef, double p, double frac, double beta) {
    return DConductivityDDepth(depth, slope, coef);
  }
};

} // namespace
} // namespace

#endif
