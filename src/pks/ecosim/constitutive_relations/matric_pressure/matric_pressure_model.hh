/*
  The hydraulic conductivity model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_MATRIC_PRESSURE_MODEL_HH_
#define AMANZI_FLOW_MATRIC_PRESSURE_MODEL_HH_

namespace Amanzi {
namespace Ecosim {
namespace Relations {

class MatricPressureModel {

 public:
  explicit MatricPressureModel(Teuchos::ParameterList& plist);

  double MatricPressure(double phi, double theta, double rho, double cv, double sl) const;

  double DMatricPressureDPorosity(double phi, double theta, double rho, double cv, double sl) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 private:

  double m_;
  double n_;
  double alpha_;
  double theta_r_;
  double sr_;
  double theta_s;

};

} //namespace
} //namespace
} //namespace

#endif
