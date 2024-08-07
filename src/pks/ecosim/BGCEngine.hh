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

#ifndef BGC_ENGINE_HH_
#define BGC_ENGINE_HH_

#include <string>
#include <vector>
#include <map>

#include "BGC_memory.hh"
#include "BGC_containers.hh"
#include "BGC_constants.hh"

#include "VerboseObject.hh"
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
  const BGCSizes& Sizes() const;

  // Initializes the data structures that hold the chemical state information.
  void InitState(BGCProperties& properties,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data,
                 int ncells_per_col_,
                 int num_components,
                 int num_columns);

  // Frees the data structures that hold the chemical state information.
  void FreeState(BGCProperties& properties,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data);
  /* Don't need for now
  void EnforceCondition(const std::string& condition_name,
                        const double time,
                        const AlquimiaProperties& properties,
                        AlquimiaState& state,
                        AlquimiaAuxiliaryData& aux_data);*/

  void DataTest();

  bool Setup(BGCProperties& properties,
               BGCState& state,
               BGCSizes& sizes,
               int num_iterations,
               int num_columns,
               int ncells_per_col_);

  // Advances the species represented by the given array of concentrations, replacing old values
  // with new values. The order of the concentrations in the array matches that of the species names
  // returned by GetSpeciesNames. Returns true if the advance is successful,
  // false if it fails.
  bool Advance(const double delta_time,
               BGCProperties& properties,
               BGCState& state,
               BGCSizes& sizes,
               int num_iterations,
               int num_columns);

  //Functions from the alquimia util section, I don't think I need the full code so I think
  //I can just copy these functions over
  void CopyBGCState(const BGCState* const source,
                         BGCState* destination);
  void CopyBGCProperties(const BGCProperties* const source,
                              BGCProperties* destination);

  //void CreateBGCInterface(const char* engine_name, BGCInterface* interface);

 private:

  // bgc data structures.
  bool bgc_initialized_;
  void* engine_state_;
  //BGCEngineFunctionality functionality_;
  BGCSizes sizes_;
  BGCInterface bgc_;
  //BGCEngineStatus bgc_status_;
  //BGCProblemMetaData bgc_metadata_;

  /*AlquimiaEngineFunctionality functionality_;
  AlquimiaSizes sizes_;
  AlquimiaInterface chem_;
  AlquimiaEngineStatus chem_status_;
  AlquimiaProblemMetaData chem_metadata_;*/

  Teuchos::RCP<VerboseObject> vo_;
  // Back-end engine name and input file.
  std::string bgc_engine_name_;
  std::string bgc_engine_inputfile_;

  //Teuchos::RCP<EcoSIM::BGCEngine> bgc_engine_;

  // forbidden.
  BGCEngine();
  BGCEngine(const BGCEngine&);
  BGCEngine& operator=(const BGCEngine&);

};

} // namespace
} // namespace

#endif
