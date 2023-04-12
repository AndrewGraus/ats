/*
  The bulk density model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ECOSIM_BULK_DENSITY_MODEL_HH_
#define AMANZI_ECOSIM_BULK_DENSITY_MODEL_HH_

namespace Amanzi {
namespace Ecosim {
namespace Relations {

class BulkDensityModel {

 public:
  explicit
  BulkDensityModel(Teuchos::ParameterList& plist);

  double BulkDensity(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;

  double DBulkDensityDPorosity(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDDensityRock(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDSaturationLiquid(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDMolarDensityLiquid(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDSaturationIce(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDMolarDensityIce(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDSaturationGas(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;
  double DBulkDensityDMolarDensityGas(double phi, double nr, double sl, double nl, double si, double ni, double sg, double ng) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif
