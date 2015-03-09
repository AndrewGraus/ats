/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_SNOW_SKIN_POTENTIAL_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SNOW_SKIN_POTENTIAL_EVALUATOR_

#include "secondary_variable_field_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class SnowSkinPotentialEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SnowSkinPotentialEvaluator(Teuchos::ParameterList& plist);

  SnowSkinPotentialEvaluator(const SnowSkinPotentialEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Key precip_key_;
  Key sd_key_;
  Key pd_key_;
  Key elev_key_;
  double factor_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SnowSkinPotentialEvaluator> factory_;

};

} //namespace
} //namespace
} //namespace

#endif