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
BulkDensityModel::BulkDensity(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return nr*(1 - phi) + phi*(ng*sg + ni*si + nl*sl);
}

double
BulkDensityModel::DBulkDensityDPorosity(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return ng*sg + ni*si + nl*sl - nr;
}

double
BulkDensityModel::DBulkDensityDDensityRock(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return 1 - phi;
}

double
BulkDensityModel::DBulkDensityDSaturationLiquid(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return nl*phi;
}

double
BulkDensityModel::DBulkDensityDMolarDensityLiquid(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return phi*sl;
}

double
BulkDensityModel::DBulkDensityDSaturationIce(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return ni*phi;
}

double
BulkDensityModel::DBulkDensityDMolarDensityIce(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return phi*si;
}

double
BulkDensityModel::DBulkDensityDSaturationGas(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return ng*phi;
}

double
BulkDensityModel::DBulkDensityDMolarDensityGas(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const
{
  return phi*sg;
}

} //namespace
} //namespace
} //namespace
