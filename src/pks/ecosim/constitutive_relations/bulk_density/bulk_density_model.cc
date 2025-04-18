/*
  The bulk density model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "bulk_density_model.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

// Constructor from ParameterList
BulkDensityModel::BulkDensityModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
BulkDensityModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{

}


// main method
double
BulkDensityModel::BulkDensity(double phi, double nr) const
{
  return nr*(1 - phi);
}

double
BulkDensityModel::DBulkDensityDPorosity(double phi, double nr) const
{
  return -1.0*nr;
}

double
BulkDensityModel::DBulkDensityDDensityRock(double phi, double nr) const
{
  return 1 - phi;
}

} //namespace
} //namespace
} //namespace
