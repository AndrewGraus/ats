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
namespace Ecosim {
namespace Relations {

// Constructor from ParameterList
MatricPressureModel::MatricPressureModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void MatricPressureModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  m_ = plist.get<double>("van genuchten m [-]", 0.2308);
  n_ = plist.get<double>("van genuchten n [-]", 1.3);
  alpha_ = plist.get<double>("van genuchten alpha [Pa^-1]", 5e-5);
  sr_ = plist.get<double>("residual saturation [-]", 0.5);
}


// main method
double MatricPressureModel::MatricPressure(double phi, double theta, double rho, double cv, double sl) const
{
  //theta_r_ = cv*rho*phi*sr_; This is water content at residual saturation
  //theta_s = cv*rho*phi; This is the water content at saturation
  //return -1.0 / std::pow(alpha_, n_) * std::pow(1.0 - std::pow((theta - cv*rho*phi*sr_) / ((cv*rho*phi) - cv*rho*phi*sr_), 1.0 / m_), -n_);
  double Se = (sl - sr_) / (1.0 - sr_);
  double Mat_p = -(1.0 / alpha_) * std::pow(std::pow(Se, -1.0 / m_) - 1.0, 1.0 / n_);
  return Mat_p;
}

double MatricPressureModel::DMatricPressureDPorosity(double phi, double theta, double rho, double cv, double sl) const
{
  return 1.0;
}

} //namespace
} //namespace
} //namespace
