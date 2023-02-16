/*
  The bulk density evaluator is an algebraic evaluator of a given model.
Richards water content evaluator: the standard form as a function of liquid saturation.  
  Generated via evaluator_generator.
*/

#include "bulk_density_evaluator.hh"
#include "bulk_density_model.hh"

namespace Amanzi {
namespace Ecosim_pk {
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
BulkDensityEvaluator::BulkDensityEvaluator(const BulkDensityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    phi_key_(other.phi_key_),
    nr_key_(other.nr_key_),
    sl_key_(other.sl_key_),
    nl_key_(other.nl_key_),
    si_key_(other.si_key_),
    ni_key_(other.ni_key_),
    sg_key_(other.sg_key_),
    ng_key_(other.ng_key_),    
    model_(other.model_) {}


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
  Key domain_name = Keys::getDomain(my_key_);

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(phi_key_);

  // dependency: density_rock
  nr_key_ = Keys::readKey(plist_, domain_name, "density rock", "density_rock");
  dependencies_.insert(nr_key_);

  // dependency: saturation_liquid
  sl_key_ = Keys::readKey(plist_, domain_name, "saturation liquid", "saturation_liquid");
  dependencies_.insert(sl_key_);

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(nl_key_);

  // dependency: saturation_ice
  si_key_ = Keys::readKey(plist_, domain_name, "saturation ice", "saturation_ice");
  dependencies_.insert(si_key_);

  // dependency: molar_density_ice
  ni_key_ = Keys::readKey(plist_, domain_name, "molar density ice", "molar_density_ice");
  dependencies_.insert(ni_key_);

  // dependency: saturation_gas
  sg_key_ = Keys::readKey(plist_, domain_name, "saturation gas", "saturation_gas");
  dependencies_.insert(sg_key_);

  // dependency: molar_density_gas
  ng_key_ = Keys::readKey(plist_, domain_name, "molar density gas", "molar_density_gas");
  dependencies_.insert(ng_key_);
}


void
BulkDensityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> nr = S->GetFieldData(nr_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> si = S->GetFieldData(si_key_);
Teuchos::RCP<const CompositeVector> ni = S->GetFieldData(ni_key_);
Teuchos::RCP<const CompositeVector> sg = S->GetFieldData(sg_key_);
Teuchos::RCP<const CompositeVector> ng = S->GetFieldData(ng_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
    const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
    const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
    const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->BulkDensity(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
    }
  }
}


void
BulkDensityEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> nr = S->GetFieldData(nr_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> si = S->GetFieldData(si_key_);
Teuchos::RCP<const CompositeVector> ni = S->GetFieldData(ni_key_);
Teuchos::RCP<const CompositeVector> sg = S->GetFieldData(sg_key_);
Teuchos::RCP<const CompositeVector> ng = S->GetFieldData(ng_key_);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDPorosity(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == nr_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDDensityRock(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDSaturationLiquid(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDMolarDensityLiquid(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == si_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDSaturationIce(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == ni_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDMolarDensityIce(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == sg_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDSaturationGas(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else if (wrt_key == ng_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& nr_v = *nr->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DBulkDensityDMolarDensityGas(phi_v[0][i], nr_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], sg_v[0][i], ng_v[0][i]);
      }
    }

  } else {
    ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace