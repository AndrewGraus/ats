/*
  The bulk density evaluator is an algebraic evaluator of a given model.
Richards water content evaluator: the standard form as a function of liquid saturation.
  Generated via evaluator_generator.
*/

#include "bulk_density_evaluator.hh"
#include "bulk_density_model.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

// Constructor from ParameterList
BulkDensityEvaluator::BulkDensityEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("bulk_density parameters");
  model_ = Teuchos::rcp(new BulkDensityModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
//Other examples don't seem to have this copy (unneccesary?)
/*BulkDensityEvaluator::BulkDensityEvaluator(const BulkDensityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    phi_key_(other.phi_key_),
    nr_key_(other.nr_key_),
    model_(other.model_) {}*/


// Virtual copy constructor
Teuchos::RCP<Evaluator>
BulkDensityEvaluator::Clone() const
{
  return Teuchos::rcp(new BulkDensityEvaluator(*this));
}


// Initialize by setting up dependencies
void
BulkDensityEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  // Seems to not handle keys or tags right fixing:
  //Key domain_name = Keys::getDomain(my_key_);
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(KeyTag{ phi_key_, tag });

  // dependency: density_rock
  nr_key_ = Keys::readKey(plist_, domain_name, "density rock", "density_rock");
  dependencies_.insert(KeyTag{ nr_key_, tag});
}


void
BulkDensityEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> nr = S.GetPtr<CompositeVector>(nr_key_, tag);

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->BulkDensity(phi_v[0][i], nr_v[0][i]);
    }
  }
}


void
BulkDensityEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> nr = S.GetPtr<CompositeVector>(nr_key_, tag);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDPorosity(phi_v[0][i], nr_v[0][i]);
      }
    }

  } else if (wrt_key == nr_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDDensityRock(phi_v[0][i], nr_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
