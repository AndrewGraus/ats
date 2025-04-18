/*
  The matric pressure evaluator is an algebraic evaluator of a given model.
Richards water content evaluator: the standard form as a function of liquid saturation.
  Generated via evaluator_generator.
*/

#include "matric_pressure_evaluator.hh"
#include "matric_pressure_model.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

// Constructor from ParameterList
MatricPressureEvaluator::MatricPressureEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("matric_pressure parameters");
  model_ = Teuchos::rcp(new MatricPressureModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
//Don't seem to need this
/*MatricPressureEvaluator::MatricPressureEvaluator(const MatricPressureEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    k_key_(other.k_key_)
    rho_key_(other.rho_key_),
    mu_key_(other.mu_key_),
    model_(other.model_) {}*/


// Virtual copy constructor
Teuchos::RCP<Evaluator>
MatricPressureEvaluator::Clone() const
{
  return Teuchos::rcp(new MatricPressureEvaluator(*this));
}


// Initialize by setting up dependencies
void
MatricPressureEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  //Key domain_name = Keys::getDomain(my_key_);
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: permeability
  porosity_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(KeyTag{ porosity_key_, tag });

  // dependency: saturation_liquid
  water_content_key_ = Keys::readKey(plist_, domain_name, "water content", "water_content");
  dependencies_.insert(KeyTag{ water_content_key_, tag});

  // dependency: molar density liquid
  mdens_liquid_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ mdens_liquid_key_, tag});

  // dependency: cell volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag});

  sl_key_ = Keys::readKey(plist_, domain_name, "liquid saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sl_key_, tag});
}


void
MatricPressureEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(porosity_key_, tag);
  Teuchos::RCP<const CompositeVector> theta = S.GetPtr<CompositeVector>(water_content_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(mdens_liquid_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& theta_v = *theta->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->MatricPressure(phi_v[0][i], theta_v[0][i], rho_v[0][i], cv_v[0][i], sl_v[0][i]);
    }
  }
}


void
MatricPressureEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(porosity_key_, tag);
  Teuchos::RCP<const CompositeVector> theta = S.GetPtr<CompositeVector>(water_content_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(mdens_liquid_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);

  if (wrt_key == porosity_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& theta_v = *theta->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);

      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMatricPressureDPorosity(phi_v[0][i], theta_v[0][i], rho_v[0][i], cv_v[0][i], sl_v[0][i]);
      }
    }

  } 
  else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
