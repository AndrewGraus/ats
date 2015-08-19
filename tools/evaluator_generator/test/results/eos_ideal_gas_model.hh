/*
  The ideal gas equation of state model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
    modelInitializeParamsList =   cv_ = plist.get<double>("heat capacity");
    evalName = eos_ideal_gas
    modelMethodDeclaration =   double Density(double temp, double pres) const;
    namespaceCaps = GENERAL
    namespace = General
    paramDeclarationList =   double cv_;
    evalNameCaps = EOS_IDEAL_GAS
    myMethodArgs = temp_v[0][i], pres_v[0][i]
    myKeyMethod = Density
    myKeyFirst = density
    evalNameString = ideal gas equation of state
    myMethodDeclarationArgs = double temp, double pres
    evalClassName = EosIdealGas
    myKey = density
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_GENERAL_EOS_IDEAL_GAS_MODEL_HH_
#define AMANZI_GENERAL_EOS_IDEAL_GAS_MODEL_HH_

namespace Amanzi {
namespace General {
namespace Relations {

class EosIdealGasModel {

 public:
  explicit
  EosIdealGasModel(Teuchos::ParameterList& plist);

  double Density(double temp, double pres) const;

  double DDensityDTemperature(double temp, double pres) const;
  double DDensityDPressure(double temp, double pres) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double cv_;

};

} //namespace
} //namespace
} //namespace

#endif