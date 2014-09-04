/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for the NullEnergy PK.  This PK simply provides a constant
   temperature, and is provided for testing with other PKs that depend upon an
   energy equation.  This could easily be provided by the state as an independent
   variable, but this is nice for testing the full hierarchy with a simple PK.

   Example usage:

   <ParameterList name="energy">
   <Parameter name="PK model" type="string" value="Constant Temperature"/>
   <Parameter name="Constant Temperature" type="double" value="290.0"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "test_snow_dist.hh"

namespace Amanzi {

void TestSnowDist::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  S->RequireField(key_, name_)->SetMesh(S->GetMesh("surface"))
          ->SetComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("precipitation_snow")->SetMesh(S->GetMesh("surface"))
      ->AddComponent("cell", AmanziMesh::CELL,1);
  S->RequireFieldEvaluator("precipitation_snow");
};

// -- call your favorite
bool TestSnowDist::advance(double dt) {

  if (sink_type_ == "factor") {
    S_next_->GetFieldData("snow_depth",name_)->Scale(sink_value_);
  } else if (sink_type_ == "constant") {
    Epetra_MultiVector& sd = *S_next_->GetFieldData("snow_depth",name_)->ViewComponent("cell",false);
    for (int c=0; c!=sd.MyLength(); ++c) {
      sd[0][c] -= 10*dt*sink_value_;
      sd[0][c] = std::max(0.,sd[0][c]);
    }
  }
    
  S_next_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetFieldData("snow_depth", name_)
    ->Update(10.*dt, *S_next_->GetFieldData("precipitation_snow"), 1.); // factor of 10 for SWE-to-snow ht conversion

  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};



} // namespace
