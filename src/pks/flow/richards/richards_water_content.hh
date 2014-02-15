/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Ethan Coon (ecoon@lanl.gov)

FieldEvaluator for water content.

Wrapping this conserved quantity as a field evaluator makes it easier to take
derivatives, keep updated, and the like.  The equation for this is simply:

WC = phi * (s_liquid * n_liquid + omega_gas * s_gas * n_gas)

This is simply the conserved quantity in Richards equation.
----------------------------------------------------------------------------- */


#ifndef AMANZI_RICHARDS_WATER_CONTENT_HH_
#define AMANZI_RICHARDS_WATER_CONTENT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class RichardsWaterContent : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  RichardsWaterContent(Teuchos::ParameterList& wc_plist);
  RichardsWaterContent(const RichardsWaterContent& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  bool is_vapor_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,RichardsWaterContent> reg_;

  
};

} // namespace
} // namespace

#endif
