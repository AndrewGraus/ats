/*
  Alquimia

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson
           Sergi Molins <smolins@lbl.gov>

  This implements the Alquimia chemistry engine.
*/

#include <iostream>
#include <cstring>
#include <cstdio>
#include <assert.h>
#include "BGCEngine.hh"
#include "errors.hh"
#include "exceptions.hh"

// Support for manipulating floating point exception handling.
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include <fenv.h>
#endif

namespace Amanzi {
namespace EcoSIM {

namespace {

//Here are the major function that the engine will need
void CopyBGCState(BGCState* dest, BGCState* src)
{
  memcpy(dest->liquid_density.data, src->liquid_density.data, sizeof(double) * src->liquid_density.size);
  memcpy(dest->gas_density.data, src->gas_density.data, sizeof(double) * src->gas_density.size);
  memcpy(dest->ice_density.data, src->ice_density.data, sizeof(double) * src->ice_density.size);
  memcpy(dest->porosity.data, src->porosity.data, sizeof(double) * src->porosity.size);
  memcpy(dest->water_content.data, src->water_content.data, sizeof(double) * src->water_content.size);
  memcpy(dest->temperature.data, src->temperature.data, sizeof(double) * src->temperature.size);
  memcpy(dest->hydraulic_conductivity.data, src->hydraulic_conductivity.data, sizeof(double) * src->hydraulic_conductivity.size);
  //memcpy(dest->total_mobile.data, src->total_mobile.data, sizeof(double) * src->total_mobile.size);

  //dest->water_density = src->water_density;
  //memcpy(dest->total_mobile.data, src->total_mobile.data, sizeof(double) * src->total_mobile.size);
}

// These functions are going into the next release of Alquimia.
void CopyBGCProperties(BGCProperties* dest, BGCProperties* src)
{
  //NEED TO EDIT STILL
  memcpy(dest->liquid_saturation.data, src->liquid_saturation.data, sizeof(double) * src->liquid_saturation.size);
  memcpy(dest->gas_saturation.data, src->gas_saturation.data, sizeof(double) * src->gas_saturation.size);
  memcpy(dest->ice_saturation.data, src->ice_saturation.data, sizeof(double) * src->ice_saturation.size);
  memcpy(dest->elevation.data, src->elevation.data, sizeof(double) * src->elevation.size);
  memcpy(dest->relative_permeability.data, src->relative_permeability.data, sizeof(double) * src->relative_permeability.size);
  memcpy(dest->thermal_conductivity.data, src->thermal_conductivity.data, sizeof(double) * src->thermal_conductivity.size);
  memcpy(dest->volume.data, src->volume.data, sizeof(double) * src->volume.size);
}

void CopyBGCAuxiliaryData(BGCAuxiliaryData* dest, BGCAuxiliaryData* src)
{
  memcpy(dest->aux_ints.data, src->aux_ints.data, sizeof(int) * src->aux_ints.size);
  memcpy(dest->aux_doubles.data, src->aux_doubles.data, sizeof(double) * src->aux_doubles.size);
}

} // namespace


BGCEngine::BGCEngine(const std::string& engineName,
                                 const std::string& inputFile) :
  bgc_engine_name_(engineName),
  bgc_engine_inputfile_(inputFile)
{
  Errors::Message msg;

  bool hands_off = false;
  if (bgc_engine_name_ == "EcoSIM") hands_off = true;

  if (bgc_engine_name_ != "EcoSIM")
  {
    std::cout << "BGCEngine: Unsupported bgc engine: '" << bgc_engine_name_ << std::endl;
    std::cout << "  Only option for now is 'EcoSIM'" << std::endl;
    Exceptions::amanzi_throw(msg);
    std::cout << "Creating the BGC Engine interface" << std::endl;
    CreateBGCInterface(bgc_engine_name_.c_str(),
                      &bgc_);

    //std::cout << "Running bgc engine setup" << std::endl;
    //bgc_.Setup(bgc_engine_inputfile_.c_str(),&sizes_);
  }

  // All alquimia function calls require a status object.
  // AllocateAlquimiaEngineStatus is in the alquimia code in alquimia_memory.c
  // The Allocate engine status function creates an object called "status"
  // and fills an element called message with calloc. This is maybe just
  // for error handling so I think it can be ignored
  //AllocateAlquimiaEngineStatus(&chem_status_);

  //CreateAlquimiaInterface - in alquimia code in alquimia_interface.c
  //This points the code to the processes for the two geochemistry codes
  // PFloTran and CrunchFlow and points the Chemistry engine functions to
  // Their corresponding functions in the drivers depending on which
  // of the geochemsitry codes. I think to start we can simply just point this
  // code directly to the EcoSIM dirver functions for the following processes:
  // Setup, Shutdown, ProcessCondition, ReactionStepOperatorSplit (advance),
  // GetAquxiliaryOutput (not neccesary for now),
  // GetProblemMetaData (not neccesary)
  // Thefore, for now I don't think we need this
  /*CreateAlquimiaInterface(chem_engine_name_.c_str(), &chem_, &chem_status_);
  if (chem_status_.error != 0)
  {
    std::cout << chem_status_.message << std::endl;
    msg << "ChemistryEngine: Could not create an interface to Alquimia.";
    Exceptions::amanzi_throw(msg);
  }*/


  //This calls the Setup function that we pointed the code to within alquimia
  /*
  chem_.Setup(chem_engine_inputfile_.c_str(),
              hands_off,
              &engine_state_,
              &sizes_,
              &functionality_,
              &chem_status_);
  if (chem_status_.error != 0)
  {
    std::cout << chem_status_.message << std::endl;
    PrintAlquimiaSizes(&sizes_, stdout);
    msg << "Error in creation of ChemistryEngine.";
    Exceptions::amanzi_throw(msg);
  }
  */
  //We don't need the chem interface call
  //Place the call to the EcoSIM setup driver here:

  /********************************
  EcoSIMSetup();
  *********************************/

  // Allocate storage for additional Alquimia data.
  // Below this sets up Alquimia metadata I don't think we
  //need this
}

BGCEngine::~BGCEngine()
{
  bgc_.Shutdown();
  //FreeAlquimiaProblemMetaData(&chem_metadata_);

  // As there are no chemical conditions, am I just deleting variables?
  /*for (GeochemicalConditionMap::iterator
       iter = chem_conditions_.begin(); iter != chem_conditions_.end(); ++iter)
  {
    //FreeBGCGeochemicalCondition(&iter->second->condition);
    FreeBGCState(&iter->second->chem_state);
    FreeBGCProperties(&iter->second->mat_props);
    FreeBGCAuxiliaryData(&iter->second->aux_data);
    delete iter->second;
  }*/

  //FreeBGCProperties(&props);
  //FreeBGCState(&state);
  //FreeBGCAuxiliaryData(&aux_data);
  //FreeAlquimiaEngineStatus(&chem_status_);
}

const BGCSizes&
BGCEngine::Sizes() const
{
  return sizes_;
}

void BGCEngine::InitState(BGCProperties& props,
                                BGCState& state,
                                BGCAuxiliaryData& aux_data,
                                int ncells_per_col_,
                                int num_components)
{
  std::cout << "Allocating Properties" << std::endl;
  std::cout << "Size of columns: " << ncells_per_col_ << std::endl;
  AllocateBGCProperties(&sizes_, &props, ncells_per_col_);
  std::cout << "Allocating State" << std::endl;
  AllocateBGCState(&sizes_, &state, ncells_per_col_, num_components);
  std::cout << "Allocating aux" << std::endl;
  AllocateBGCAuxiliaryData(&sizes_, &aux_data, ncells_per_col_);
  //AllocateAlquimiaAuxiliaryOutputData(&sizes_, &aux_output);

  // Make sure the auxiliary ints/doubles are zeroed out.
  std::fill(aux_data.aux_ints.data, aux_data.aux_ints.data + aux_data.aux_ints.size, 0);
  std::fill(aux_data.aux_doubles.data, aux_data.aux_doubles.data + aux_data.aux_doubles.size, 0.0);
}

void BGCEngine::FreeState(BGCProperties& props,
                                BGCState& state,
                                BGCAuxiliaryData& aux_data)
{
  FreeBGCProperties(&props);
  FreeBGCState(&state);
  FreeBGCAuxiliaryData(&aux_data);
  //FreeAlquimiaAuxiliaryOutputData(&aux_output);
}

/*void BGCEngine::DataTest() {

  std::cout << "Data test for calling function" << std::endl;
  bgc_.DataTest();
}*/

bool BGCEngine::Setup(BGCProperties& props,
                              BGCState& state,
                              BGCAuxiliaryData& aux_data,
                              int num_iterations,
                              int ncol)
{
  std::cout << "Running BGC Engine Setup" << std::endl;
  // Advance the chemical reaction all operator-split-like.
  bgc_.Setup(&props,
                &state,
                &aux_data,
                num_iterations,
                ncol);

  //This is alquimia's advance function which we won't need
  //calling EcoSIM advance driver

  /*************************
  EcoSIMAdvance();
  *************************/
}

bool BGCEngine::Advance(const double delta_time,
                              BGCProperties& props,
                              BGCState& state,
                              BGCAuxiliaryData& aux_data,
                              int num_iterations,
                              int ncol)
{
  std::cout << "Running BGC Engine Advance" << std::endl;
  // Advance the chemical reaction all operator-split-like.
  bgc_.Advance(delta_time,
                &props,
                &state,
                &aux_data,
                num_iterations,
                ncol);

  //This is alquimia's advance function which we won't need
  //calling EcoSIM advance driver

  /*************************
  EcoSIMAdvance();
  *************************/
}


/*void CreateBGCInterface(const char* const engine_name, BGCInterface* interface)
 {

   interface->Setup = &ecosim_setup;
   interface->Shutdown = &ecosim_shutdown;
   interface->Advance = &ecosim_advance;

 }
*/
//For now I don't need any of the rest of this code yet. Just commenting
//out in case I need it later.

/*
void ChemistryEngine::CreateCondition(const std::string& condition_name)
{
  // NOTE: a condition with zero aqueous/mineral constraints is assumed to be defined in
  // NOTE: the backend engine's input file.
  GeochemicalConditionData* condition = new GeochemicalConditionData();
  condition->processed = false;
  int num_aq = 0, num_min = 0;
  AllocateAlquimiaProperties(&sizes_, &condition->mat_props);
  AllocateAlquimiaGeochemicalCondition(kAlquimiaMaxStringLength, num_aq, num_min, &condition->condition);
  AllocateAlquimiaState(&sizes_, &condition->chem_state);
  AllocateAlquimiaAuxiliaryData(&sizes_, &condition->aux_data);
  std::strcpy(condition->condition.name, condition_name.c_str());

  // Add this to the conditions map.
  chem_conditions_[condition_name] = condition;
}

void ChemistryEngine::AddAqueousConstraint(const std::string& condition_name,
                                           const std::string& primary_species_name,
                                           const std::string& constraint_type,
                                           const std::string& associated_species)
{
  assert(condition_name.length() > 0);
  assert(primary_species_name.length() > 0);
  assert((constraint_type == "total_aqueous") || (constraint_type == "charge") ||
         (constraint_type == "free") || (constraint_type == "mineral") ||
         (constraint_type == "gas") || (constraint_type == "pH"));

  GeochemicalConditionMap::iterator iter = chem_conditions_.find(condition_name);
  if (iter != chem_conditions_.end())
  {
    AlquimiaGeochemicalCondition* condition = &iter->second->condition;

    // Do we have an existing constraint?
    int index = -1;
    for (int i = 0; i < condition->aqueous_constraints.size; ++i)
    {
      if (!std::strcmp(condition->aqueous_constraints.data[i].primary_species_name, primary_species_name.c_str()))
      {
        // Overwrite the old constraint.
        index = i;
        free(condition->aqueous_constraints.data[index].primary_species_name);
        free(condition->aqueous_constraints.data[index].constraint_type);
        if (condition->aqueous_constraints.data[index].associated_species != NULL)
          free(condition->aqueous_constraints.data[index].associated_species);
      }
    }
    if (index == -1)
    {
      // New constraint!
      index = condition->aqueous_constraints.size;
      condition->aqueous_constraints.size++;
      condition->aqueous_constraints.data = (AlquimiaAqueousConstraint*)realloc(condition->aqueous_constraints.data, sizeof(AlquimiaAqueousConstraint) * (index+1));
    }

    // Add the aqueous constraint.
    condition->aqueous_constraints.data[index].primary_species_name = strdup(primary_species_name.c_str());
    condition->aqueous_constraints.data[index].constraint_type = strdup(constraint_type.c_str());
    if (!associated_species.empty())
      condition->aqueous_constraints.data[index].associated_species = strdup(associated_species.c_str());
    else
      condition->aqueous_constraints.data[index].associated_species = NULL;
  }
  else
  {
    Errors::Message msg;
    msg << "ChemistryEngine::AddAqueousConstraint: no condition named '" << condition_name << "'.";
    Exceptions::amanzi_throw(msg);
  }
}

void ChemistryEngine::EnforceCondition(const std::string& condition_name,
                                       const double time,
                                       const AlquimiaProperties& mat_props,
                                       AlquimiaState& chem_state,
                                       AlquimiaAuxiliaryData& aux_data,
                                       AlquimiaAuxiliaryOutputData& aux_output)
{
  Errors::Message msg;

  // Retrieve the chemical condition for the given name.
  GeochemicalConditionMap::iterator iter = chem_conditions_.find(condition_name);
  if (iter == chem_conditions_.end())
  {
    CreateCondition(condition_name);
    iter = chem_conditions_.find(condition_name);
  }

#ifdef AMANZI_USE_FENV
  // Disable divide-by-zero floating point exceptions.
  int fpe_mask = fedisableexcept(FE_DIVBYZERO);
#endif

  AlquimiaGeochemicalCondition* condition = &iter->second->condition;
  AlquimiaProperties& nc_mat_props = const_cast<AlquimiaProperties&>(mat_props);
  if (!iter->second->processed)
  {
    // Copy the given state data into place for this condition.
    CopyAlquimiaProperties(&iter->second->mat_props, &nc_mat_props);
    CopyAlquimiaState(&iter->second->chem_state, &chem_state);
    CopyAlquimiaAuxiliaryData(&iter->second->aux_data, &aux_data);

    // Process the condition on the given array at the given time.
    // FIXME: Time is ignored for the moment.
    chem_.ProcessCondition(&engine_state_, condition, &iter->second->mat_props,
                           &iter->second->chem_state, &iter->second->aux_data, &chem_status_);
    iter->second->processed = true;
  }

  // Copy the constraint's data into place.
  // Here, condition names are used but Alquimia properties are
  // given as "Material" zones in Amanzi's inputthus disconnecting them
  // from "Initial" conditions (i.e. condition names)
  // For the sake of performance, we avoid here to re-process conditions
  // but we need to retain materical properties as provided in Amanzi inputs
  // CopyAlquimiaProperties(&nc_mat_props, &iter->second->mat_props);
  CopyAlquimiaState(&chem_state, &iter->second->chem_state);
  CopyAlquimiaAuxiliaryData(&aux_data, &iter->second->aux_data);

#ifdef AMANZI_USE_FENV
  // Re-enable pre-existing floating point exceptions.
  feclearexcept(fpe_mask);
  fpe_mask = feenableexcept(fpe_mask);
#endif

// FIXME: Figure out a neutral parallel-friendly way to report errors.
  assert(chem_status_.error == 0);

#if 0
  if (chem_status_.error != 0)
    ierr = -1;

  // figure out if any of the processes threw an error, if so all processes will re-throw
  int recv = 0;
  mesh_->get_comm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0)
  {
    msg << "Error in enforcement of chemical condition '" << condition_name << "'";
    Exceptions::amanzi_throw(msg);
  }
#endif
}

const std::string& ChemistryEngine::Name() const
{
  return chem_engine_name_;
}

bool ChemistryEngine::IsThreadSafe() const
{
  Errors::Message msg;

  if (not chem_initialized_)
  {
    msg << "ChemistryEngine: Cannot query before initialization!";
    Exceptions::amanzi_throw(msg);
  }

  return functionality_.thread_safe;
}

int ChemistryEngine::NumPrimarySpecies() const
{
  return sizes_.num_primary;
}

int ChemistryEngine::NumAqueousComplexes() const
{
  return sizes_.num_aqueous_complexes;
}

int ChemistryEngine::NumSorbedSpecies() const
{
  return sizes_.num_sorbed;
}

void ChemistryEngine::GetPrimarySpeciesNames(std::vector<std::string>& species_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->primary_names.size;
  species_names.resize(N);
  for (int i = 0; i < N; ++i)
    species_names[i] = std::string(metadata->primary_names.data[i]);
}

int ChemistryEngine::NumMinerals() const
{
  return sizes_.num_minerals;
}

void ChemistryEngine::GetMineralNames(std::vector<std::string>& mineral_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->mineral_names.size;
  mineral_names.resize(N);
  for (int i = 0; i < N; ++i)
    mineral_names[i] = std::string(metadata->mineral_names.data[i]);
}

int ChemistryEngine::NumSurfaceSites() const
{
  return sizes_.num_surface_sites;
}

void ChemistryEngine::GetSurfaceSiteNames(std::vector<std::string>& site_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->surface_site_names.size;
  site_names.resize(N);
  for (int i = 0; i < N; ++i)
    site_names[i] = std::string(metadata->surface_site_names.data[i]);
}

int ChemistryEngine::NumIonExchangeSites() const
{
  return sizes_.num_ion_exchange_sites;
}

void ChemistryEngine::GetIonExchangeNames(std::vector<std::string>& ion_exchange_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->ion_exchange_names.size;
  ion_exchange_names.resize(N);
  for (int i = 0; i < N; ++i)
    ion_exchange_names[i] = std::string(metadata->ion_exchange_names.data[i]);
}

int ChemistryEngine::NumIsothermSpecies() const
{
  return sizes_.num_isotherm_species;
}

void ChemistryEngine::GetIsothermSpeciesNames(std::vector<std::string>& species_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->isotherm_species_names.size;
  species_names.resize(N);
  for (int i = 0; i < N; ++i)
    species_names[i] = std::string(metadata->isotherm_species_names.data[i]);
}

int ChemistryEngine::NumFreeIonSpecies() const
{
  return sizes_.num_primary;
}

void ChemistryEngine::GetAuxiliaryOutputNames(std::vector<std::string>& aux_names,
        std::vector<std::vector<std::string>>& subfield_names) const
{
  aux_names.clear();
  subfield_names.clear();
  aux_names.emplace_back("pH");
  subfield_names.emplace_back(std::vector<std::string>{"0"});

  const AlquimiaProblemMetaData* metadata = &chem_metadata_;

  // Mineral data -- one per mineral.
  int N = metadata->mineral_names.size;
  if (N > 0) {
    std::vector<std::string> mineral_names;
    for (int i = 0; i < N; ++i) {
      mineral_names.emplace_back(metadata->mineral_names.data[i]);
    }
    aux_names.emplace_back("mineral_saturation_index");
    subfield_names.emplace_back(mineral_names);
    aux_names.emplace_back("mineral_reaction_rate");
    subfield_names.emplace_back(std::move(mineral_names));
  }

  // Auxiliary data per primary species.
  N = metadata->primary_names.size;
  if (N > 0) {
    std::vector<std::string> primary_names;
    for (int i = 0; i < N; ++i) {
      primary_names.emplace_back(metadata->primary_names.data[i]);
    }
    aux_names.emplace_back("primary_free_ion_concentration");
    subfield_names.emplace_back(primary_names);
    aux_names.emplace_back("primary_activity_coeff");
    subfield_names.emplace_back(std::move(primary_names));
  }

  // Secondary auxiliary data.
  N = this->NumAqueousComplexes();
  if (N > 0) {
    std::vector<std::string> secondary_names;
    for (int i = 0; i < N; ++i) {
      secondary_names.emplace_back(std::to_string(i));
    }
    aux_names.emplace_back("secondary_free_ion_concentration");
    subfield_names.emplace_back(secondary_names);
    aux_names.emplace_back("secondary_activity_coeff");
    subfield_names.emplace_back(std::move(secondary_names));
  }
}

int ChemistryEngine::NumAqueousKinetics() const
{
  return sizes_.num_aqueous_kinetics;
}

void ChemistryEngine::GetAqueousKineticNames(std::vector<std::string>& kinetics_names) const
{
  const AlquimiaProblemMetaData* metadata = &chem_metadata_;
  int N = metadata->aqueous_kinetic_names.size;
  kinetics_names.resize(N);
  for (int i = 0; i < N; ++i)
    kinetics_names[i] = std::string(metadata->aqueous_kinetic_names.data[i]);
}

const AlquimiaSizes& ChemistryEngine::Sizes() const
{
  return sizes_;
}*/

} // namespace
} // namespace
