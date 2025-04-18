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

  double BulkDensity(double phi, double nr) const;

  double DBulkDensityDPorosity(double phi, double nr) const;
  double DBulkDensityDDensityRock(double phi, double nr) const;

 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif
