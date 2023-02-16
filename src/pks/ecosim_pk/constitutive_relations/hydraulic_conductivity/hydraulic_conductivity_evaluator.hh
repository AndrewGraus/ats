/*
  The hydraulic conductivity evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ECOSIM_PK_HYDRAULIC_CONDUCTIVITY_EVALUATOR_HH_
#define AMANZI_ECOSIM_PK_HYDRAULIC_CONDUCTIVITY_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Ecosim_pk {
namespace Relations {

class HydraulicConductivityModel;

class HydraulicConductivityEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  HydraulicConductivityEvaluator(Teuchos::ParameterList& plist);
  HydraulicConductivityEvaluator(const HydraulicConductivityEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<HydraulicConductivityModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key k_key_;
  Key rho_key_;
  Key mu_key_;

  Teuchos::RCP<HydraulicConductivityModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,HydraulicConductivityEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif