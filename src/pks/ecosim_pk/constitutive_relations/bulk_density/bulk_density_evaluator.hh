/*
  The bulk density evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ECOSIM_PK_BULK_DENSITY_EVALUATOR_HH_
#define AMANZI_ECOSIM_PK_BULK_DENSITY_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Ecosim_pk {
namespace Relations {

class BulkDensityModel;

class BulkDensityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit
  BulkDensityEvaluator(Teuchos::ParameterList& plist);
  BulkDensityEvaluator(const BulkDensityEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<BulkDensityModel> get_model() { return model_; }

 protected:
   // Required methods from EvaluatorSecondaryMonotypeCV
   virtual void Evaluate_(const State& S,
           const std::vector<CompositeVector*>& result) override;
   virtual void EvaluatePartialDerivative_(const State& S,
           const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

   void InitializeFromPlist_();

 protected:
  Key phi_key_;
  Key nr_key_;
  Key sl_key_;
  Key nl_key_;
  Key si_key_;
  Key ni_key_;
  Key sg_key_;
  Key ng_key_;

  Teuchos::RCP<BulkDensityModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,BulkDensityEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
