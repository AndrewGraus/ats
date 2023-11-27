/*
  The hydraulic conductivity evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_MATRIC_PRESSURE_EVALUATOR_HH_
#define AMANZI_FLOW_MATRIC_PRESSURE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

class MatricPressureModel;

class MatricPressureEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit MatricPressureEvaluator(Teuchos::ParameterList& plist);
  MatricPressureEvaluator(const MatricPressureEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<MatricPressureModel> get_model() { return model_; }

 protected:
   // Required methods from EvaluatorSecondaryMonotypeCV
   virtual void Evaluate_(const State& S,
           const std::vector<CompositeVector*>& result) override;
   virtual void EvaluatePartialDerivative_(const State& S,
           const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key porosity_key_;
  Key water_content_key_;
  Key mdens_liquid_key_;
  Key cv_key_;

  Teuchos::RCP<MatricPressureModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,MatricPressureEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
