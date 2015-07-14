/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for WRM implementations.
   ------------------------------------------------------------------------- */

#include "wrm_factory.hh"

// explicity instantitate the static data of Factory<WRM>
template<> 
Amanzi::Utils::Factory<Amanzi::Flow::FlowRelations::WRM>::map_type* 
Amanzi::Utils::Factory<Amanzi::Flow::FlowRelations::WRM>::map_;

