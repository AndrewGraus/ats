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
  BulkDensityEvaluator(const BulkDensityEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<BulkDensityModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

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