/*
  The hydraulic conductivity model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "matric_pressure_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
MatricPressureModel::MatricPressureModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void MatricPressureModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  m_ = plist.get<double>("van genuchten m");
  n_ = plist.get<double>("van genuchten n");
  alpha_ = plist.get<double>("van genuchten alpha");
  sr_ = plist.get<double>("residual saturation [-]");
}


// main method
double MatricPressureModel::MatricPressure(double phi, double theta, double rho, double cv) const
{
  //theta_r_ = cv*rho*phi*sr_; This is water content at residual saturation
  //theta_s = cv*rho*phi; This is the water content at saturation
  return -1.0 / std::pow(alpha_, n_) * std::pow(1.0 - std::pow((theta - cv*rho*phi*sr_) / ((cv*rho*phi) - cv*rho*phi*sr_), 1.0 / m_), -n_);
}

double MatricPressureModel::DMatricPressureDPorosity(double phi, double theta, double rho, double cv) const
{
  return 1.0;
}

} //namespace
} //namespace
} //namespace