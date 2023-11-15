/*
  The hydraulic conductivity model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_MATRIC_PRESSURE_MODEL_HH_
#define AMANZI_FLOW_MATRIC_PRESSURE_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class MatricPressureModel {

 public:
  explicit
  MatricPressureModel(Teuchos::ParameterList& plist);

  double MatricPressure(double phi, double theta, double rho, double cv) const;

  double DMatricPressureDPorosity(double phi, double theta, double rho, double cv) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double m_, n_, alpha_, theta_r_;

};

} //namespace
} //namespace
} //namespace

#endif
