/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_THREEPHASE_FACTORY_HH_
#define PK_ENERGY_RELATIONS_TC_THREEPHASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "thermal_conductivity_threephase.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class ThermalConductivityThreePhaseFactory : public Utils::Factory<ThermalConductivityThreePhase> {

public:
  Teuchos::RCP<ThermalConductivityThreePhase> createThermalConductivityModel(Teuchos::ParameterList& plist);
};

} // namespace
} // namespace
} // namespace

#endif
