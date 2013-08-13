/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_SCHEME_
#define AMANZI_UPWINDING_SCHEME_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

// forward declaration
class State;

namespace Operators {

enum UpwindMethod {
  UPWIND_METHOD_CENTERED = 0,
  UPWIND_METHOD_GRAVITY,
  UPWIND_METHOD_TOTAL_FLUX,
  UPWIND_METHOD_ARITHMETIC_MEAN,
  UPWIND_METHOD_POTENTIAL_DIFFERENCE
};

class Upwinding {

 public:
  virtual void
  Update(const Teuchos::Ptr<State>& S) = 0;
};

} // namespace
} // namespace

#endif
