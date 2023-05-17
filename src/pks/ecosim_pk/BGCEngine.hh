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
/*#include "alquimia/alquimia_util.h"
#include "alquimia/alquimia_constants.h"
#include "alquimia/alquimia_containers.h"
#include "alquimia/alquimia_interface.h"*/

namespace Amanzi {
namespace EcoSIM {

class BGCEngine {
 public:

   //This is the main struct which will hold the data for the engine including
   //its set up and shutdown regions
   typedef struct {
     /* read data files/structures, initialize memory, basis management
        (includes reading database, swapping basis, etc.) */
     void (*Setup)(
         const char* input_filename,
         //bool hands_off,
         //void* pft_engine_state,
         BGCSizes* sizes);

     /* gracefully shutdown the engine, cleanup memory */
     void (*Shutdown)(
       void* pft_engine_state);

     /* constrain processing for boundary/initial constraints. Called
        once for each IC/BC. */
     /*void (*Setup)(
         void* pft_engine_state,
         BGCGeochemicalCondition* condition,
         BGCProperties* props,
         BGCState* state,
         BGCAuxiliaryData* aux_data,
         BGCEngineStatus* status);*/

     /* take one (or more?) reaction steps in operator split mode */
     void (*Advance)(
         void* pft_engine_state,
         double delta_t,
         BGCProperties* props,
         BGCState* state,
         BGCAuxiliaryData* aux_data);

     /* Access to user selected geochemical data for output, i.e. pH,
        mineral SI, reaction rates */
     /*void (*GetAuxiliaryOutput)(
         void* pft_engine_state,
         BGCProperties* props,
         BGCState* state,
         BGCAuxiliaryData* aux_data,
         BGCAuxiliaryOutputData* aux_out,
         BGCEngineStatus* status);*/
   } BGCInterface;

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
  void InitState(BGCProperties& props,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data,
                 int ncells_per_col_,
                 int num_components);

  // Frees the data structures that hold the chemical state information.
  void FreeState(BGCProperties& props,
                 BGCState& state,
                 BGCAuxiliaryData& aux_data);
  /* Don't need for now
  void EnforceCondition(const std::string& condition_name,
                        const double time,
                        const AlquimiaProperties& props,
                        AlquimiaState& state,
                        AlquimiaAuxiliaryData& aux_data);*/

  // Advances the species represented by the given array of concentrations, replacing old values
  // with new values. The order of the concentrations in the array matches that of the species names
  // returned by GetSpeciesNames. Returns true if the advance is successful,
  // false if it fails.
  bool Advance(const double delta_time,
               const BGCProperties& props,
               BGCState& state,
               BGCAuxiliaryData& aux_data,
               int& num_iterations);

  //Functions from the alquimia util section, I don't think I need the full code so I think
  //I can just copy these functions over
  void CopyBGCState(const BGCState* const source,
                         BGCState* destination);
  void CopyBGCProperties(const BGCProperties* const source,
                              BGCProperties* destination);

  void CreateBGCInterface(const char* const engine_name,
                          BGCInterface* interface);

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
