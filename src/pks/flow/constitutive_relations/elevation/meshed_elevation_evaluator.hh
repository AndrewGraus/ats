/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MESHED_ELEVATION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_MESHED_ELEVATION_EVALUATOR_

#include "factory.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class MeshedElevationEvaluator : public ElevationEvaluator {

 public:
  explicit
  MeshedElevationEvaluator(Teuchos::ParameterList& plist);

  MeshedElevationEvaluator(const MeshedElevationEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 private:
  static Utils::RegisteredFactory<FieldEvaluator,MeshedElevationEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
