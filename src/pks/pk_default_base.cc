/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "pk_default_base.hh"

namespace Amanzi {


void PKDefaultBase::setup(const Teuchos::Ptr<State>& S) {
  // THIS MAY BE CALLED MORE THAN ONCE!
  name_ = plist_.get<std::string>("PK name");

  // set up the VerboseObject
  vo_ = Teuchos::rcp(new VerboseObject(name_, plist_));
}


void PKDefaultBase::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


void PKDefaultBase::solution_to_state(const Teuchos::RCP<const TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  Teuchos::RCP<TreeVector> nc_soln = Teuchos::rcp_const_cast<TreeVector>(soln);
  solution_to_state(nc_soln, S);
}

} // namespace
