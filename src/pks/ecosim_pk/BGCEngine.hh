/*
  Alqumia

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson

  This is a point of contact for the chemistry engine exposed by Alquimia
  to the rest of Amanzi--it provides the ability to enforce geochemical
  conditions and to integrate reactions given a chemical configuration.
*/

#ifndef AMANZI_CHEMISTRY_ENGINE_HH_
#define AMANZI_CHEMISTRY_ENGINE_HH_

#include <string>
#include <vector>
#include <map>

#include "BGC_memory.h"
/*#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_containers.h"
#include "alquimia/alquimia_interface.h"*/

namespace Amanzi {
namespace EcoSIM {

class BGCEngine {

 public:

  // Constructs a chemistry engine using the given engine (backend) name and input file.
  BGCEngine(const std::string& engineName, const std::string& inputFile);

  // Destructor.
  ~BGCEngine();

  // Returns the name of the backend that does the chemistry.
  const std::string& Name() const;

  // Returns true if the chemistry engine is thread-safe, false if not.
  bool IsThreadSafe() const;

  // Returns a reference to a "sizes" object that can be queried to find the sizes of the various
  // arrays representing the geochemical state within the engine.
  const AlquimiaSizes& Sizes() const;

  // Initializes the data structures that hold the chemical state information.
  void InitState(BGCProperties& mat_props,
                 BGCState& chem_state,
                 BGCAuxiliaryData& aux_data);

  // Frees the data structures that hold the chemical state information.
  void FreeState(BGCProperties& mat_props,
                 BGCState& chem_state,
                 BGCAuxiliaryData& aux_data);

  void EnforceCondition(const std::string& condition_name,
                        const double time,
                        const AlquimiaProperties& mat_props,
                        AlquimiaState& chem_state,
                        AlquimiaAuxiliaryData& aux_data,
                        AlquimiaAuxiliaryOutputData& aux_output);

  // Advances the species represented by the given array of concentrations, replacing old values
  // with new values. The order of the concentrations in the array matches that of the species names
  // returned by GetSpeciesNames. Returns true if the advance is successful,
  // false if it fails.
  bool Advance(const double delta_time,
               const BGCProperties& mat_props,
               BGCState& chem_state,
               BGCAuxiliaryData& aux_data,
               int& num_iterations);

 private:

  // bgc data structures.
  bool bgc_initialized_;
  void* engine_state_;

  /*AlquimiaEngineFunctionality functionality_;
  AlquimiaSizes sizes_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaProblemMetaData chem_metadata_;*/

  // Back-end engine name and input file.
  std::string bgc_engine_name_;
  std::string bgc_engine_inputfile_;

  // forbidden.
  BGCEngine();
  BGCEngine(const BGCEngine&);
  BGCEngine& operator=(const BGCEngine&);

};

} // namespace
} // namespace

#endif
