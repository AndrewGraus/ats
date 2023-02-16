/*
  The hydraulic conductivity model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "hydraulic_conductivity_model.hh"

namespace Amanzi {
namespace Ecosim_pk {
namespace Relations {

// Constructor from ParameterList
HydraulicConductivityModel::HydraulicConductivityModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
HydraulicConductivityModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  g_ = plist.get<double>("gravitational constant");
}


// main method
double
HydraulicConductivityModel::HydraulicConductivity(double k, double rho, double mu) const
{
  return k*rho;
}

double
HydraulicConductivityModel::DHydraulicConductivityDPermeabiilty(double k, double rho, double mu) const
{
  return rho;
}

double
HydraulicConductivityModel::DHydraulicConductivityDMassDensityLiquid(double k, double rho, double mu) const
{
  return k;
}

double
HydraulicConductivityModel::DHydraulicConductivityDViscosityLiquid(double k, double rho, double mu) const
{
  return 0;
}

} //namespace
} //namespace
} //namespace
  