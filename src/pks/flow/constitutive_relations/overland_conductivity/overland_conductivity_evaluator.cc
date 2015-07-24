/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "independent_variable_field_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "manning_conductivity_model.hh"
#include "split_denominator_conductivity_model.hh"
#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

OverlandConductivityEvaluator::OverlandConductivityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  depth_key_ = plist_.get<std::string>("height key", "ponded_depth");
  dependencies_.insert(depth_key_);

  slope_key_ = plist_.get<std::string>("slope key", "slope_magnitude");
  dependencies_.insert(slope_key_);

  coef_key_ = plist_.get<std::string>("coefficient key", "manning_coefficient");
  dependencies_.insert(coef_key_);

  dt_ = plist_.get<bool>("include dt factor", false);
  if (dt_) {
    factor_ = plist_.get<double>("dt factor");
  }
  dens_ = plist_.get<bool>("include density factor", true);

  if (dens_) {
    dens_key_ = plist_.get<std::string>("density key", "surface-molar_density_liquid");
    dependencies_.insert(dens_key_);
  }
  
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("overland conductivity key",
            "overland_conductivity");
  }

  // create the model
  ASSERT(plist_.isSublist("overland conductivity model"));
  Teuchos::ParameterList sublist = plist_.sublist("overland conductivity model");
  std::string model_type = sublist.get<std::string>("overland conductivity type", "manning");
  if ((model_type == "manning") || (model_type == "harmonic mean manning")){
    model_ = Teuchos::rcp(new ManningConductivityModel(sublist));
  } else if (model_type == "split denominator") {
    model_ = Teuchos::rcp(new SplitDenominatorConductivityModel(sublist));
  } else if (model_type == "ponded depth passthrough") {
    model_ = Teuchos::rcp(new PondedDepthPassthroughConductivityModel(sublist));
  } else {
    ASSERT(0);
  }
    
}


OverlandConductivityEvaluator::OverlandConductivityEvaluator(const OverlandConductivityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    depth_key_(other.depth_key_),
    slope_key_(other.slope_key_),
    coef_key_(other.coef_key_),
    dens_key_(other.dens_key_),
    dt_(other.dt_),
    factor_(other.factor_),
    dens_(other.dens_),
    model_(other.model_) {}


Teuchos::RCP<FieldEvaluator>
OverlandConductivityEvaluator::Clone() const {
  return Teuchos::rcp(new OverlandConductivityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void OverlandConductivityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
    const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
    const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    if (dt_) {
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      }
    } else {
      //      double dt = *S->GetScalarData("dt");
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
      }
    }

    if (dens_) {
      const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(*comp,false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] *= dens_v[0][i];
      }
    }
  }
}


void OverlandConductivityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetFieldData(depth_key_);
  Teuchos::RCP<const CompositeVector> slope = S->GetFieldData(slope_key_);
  Teuchos::RCP<const CompositeVector> coef = S->GetFieldData(coef_key_);

  if (wrt_key == depth_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      if (dt_) {
        //        double dt = *S->GetScalarData("dt");
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]) * factor_;
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->DConductivityDDepth(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }

      if (dens_) {
        const Epetra_MultiVector& dens_v = *S->GetFieldData(dens_key_)->ViewComponent(*comp,false);
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] *= dens_v[0][i];
        }
      }
    }

  } else if (wrt_key == dens_key_) {
    ASSERT(dens_);
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& depth_v = *depth->ViewComponent(*comp,false);
      const Epetra_MultiVector& slope_v = *slope->ViewComponent(*comp,false);
      const Epetra_MultiVector& coef_v = *coef->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      if (dt_) {
        //        double dt = *S->GetScalarData("dt");
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(factor_*depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      } else {
        for (int i=0; i!=ncomp; ++i) {
          result_v[0][i] = model_->Conductivity(depth_v[0][i], slope_v[0][i], coef_v[0][i]);
        }
      }
    }

  } else {
    ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace

