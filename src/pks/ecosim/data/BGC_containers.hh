/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
** Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
** through Lawrence Berkeley National Laboratory (subject to receipt of any
** required approvals from the U.S. Dept. of Energy).  All rights reserved.
**
** Alquimia is available under a BSD license. See LICENSE.txt for more
** information.
**
** If you have questions about your rights to use or distribute this software,
** please contact Berkeley Lab's Technology Transfer and Intellectual Property
** Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
**
** NOTICE.  This software was developed under funding from the U.S. Department
** of Energy.  As such, the U.S. Government has been granted for itself and
** others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
** license in the Software to reproduce, prepare derivative works, and perform
** publicly and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of Energy,
** and subject to any subsequent five (5) year renewals, the U.S. Government is
** granted for itself and others acting on its behalf a paid-up, nonexclusive,
** irrevocable, worldwide license in the Software to reproduce, prepare derivative
** works, distribute copies to the public, perform publicly and display publicly,
** and to permit others to do so.
**
** Authors: Benjamin Andre <bandre@lbl.gov>
*/

#ifndef ALQUIMIA_CONTAINERS_H_
#define ALQUIMIA_CONTAINERS_H_

/*******************************************************************************
 **
 ** C implementation of the alquimia containers.
 **
 ** These are passed directly into the fortran routines. The
 ** signatures must match exactly with the fortran side of things.
 **
 ******************************************************************************/

//#include "alquimia/alquimia.h"
//#include <cstddef>
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  typedef struct {
    int size, capacity;
    double* data;
  } BGCVectorDouble;

  typedef struct {
    int size, capacity;
    int* data;
  } BGCVectorInt;

  typedef struct {
    /* NOTE: this is a vector of strings */
    int size, capacity;
    char** data;
  } BGCVectorString;

  typedef struct {
    int cells, columns, capacity_cells, capacity_columns;
    double* data;
  } BGCMatrixDouble;

  typedef struct {
    int cells, columns, capacity_cells, capacity_columns;
    int* data;
  } BGCMatrixInt;

  typedef struct {
    int cells, columns, capacity;
    char* data;
  } BGCMatrixString;

  typedef struct {
    int cells, columns, components, capacity_cells, capacity_columns, capacity_components;
    double*** data;
  } BGCTensorDouble;

  typedef struct {
    int cells, columns, components, capacity_cells, capacity_columns, capacity_components;
    int*** data;
  } BGCTensorInt;

  typedef struct {
    int cells, columns, components, capacity;
    char*** data;
  } BGCTensorString;

  typedef struct {
    int ncells_per_col_;
    int num_components;
    int num_columns;
  } BGCSizes;

  typedef struct {
    BGCMatrixDouble liquid_density;
    BGCMatrixDouble gas_density;
    BGCMatrixDouble ice_density;
    BGCMatrixDouble porosity;
    BGCMatrixDouble water_content;
    BGCMatrixDouble matric_pressure;
    BGCMatrixDouble temperature;
    BGCMatrixDouble hydraulic_conductivity;
    BGCMatrixDouble bulk_density;
    BGCMatrixDouble subsurface_water_source;
    BGCMatrixDouble subsurface_energy_source;
    BGCVectorDouble surface_energy_source;
    BGCVectorDouble surface_water_source;
    BGCVectorDouble snow_depth;
    BGCVectorDouble canopy_height;
    BGCTensorDouble total_component_concentration;
  } BGCState;

  typedef struct {
    BGCMatrixDouble liquid_saturation;
    BGCMatrixDouble gas_saturation;
    BGCMatrixDouble ice_saturation;
    BGCMatrixDouble relative_permeability;
    BGCMatrixDouble thermal_conductivity;
    BGCMatrixDouble volume;
    BGCMatrixDouble depth;
    BGCMatrixDouble depth_c;
    BGCMatrixDouble dz;
    BGCMatrixDouble plant_wilting_factor;
    BGCMatrixDouble rooting_depth_fraction;
    BGCVectorDouble shortwave_radiation;
    BGCVectorDouble longwave_radiation;
    BGCVectorDouble air_temperature;
    BGCVectorDouble vapor_pressure_air;
    BGCVectorDouble wind_speed;
    BGCVectorDouble precipitation;
    BGCVectorDouble precipitation_snow;
    BGCVectorDouble elevation;
    BGCVectorDouble aspect;
    BGCVectorDouble slope;
    BGCVectorDouble LAI;
    BGCVectorDouble SAI;
    BGCVectorDouble vegetation_type;
    double atm_n2;
    double atm_o2;
    double atm_co2;
    double atm_ch4;
    double atm_n2o;
    double atm_h2;
    double atm_nh3;
    double heat_capacity;
    double field_capacity;
    double wilting_point;
    int current_day;
    int current_year;
    bool p_bool;
  } BGCProperties;

  /*
  typedef struct {
    BGCVectorDouble liquid_density;
    BGCVectorDouble gas_density;
    BGCVectorDouble ice_density;
    BGCVectorDouble porosity;
    BGCVectorDouble water_content;
    BGCVectorDouble suction_head;
    BGCVectorDouble temperature;
    BGCVectorDouble hydraulic_conductivity;
    BGCVectorDouble bulk_density;
    BGCMatrixDouble total_component_concentration;
  } BGCState;

  typedef struct {
    BGCVectorDouble liquid_saturation;
    BGCVectorDouble gas_saturation;
    BGCVectorDouble ice_saturation;
    BGCVectorDouble relative_permeability;
    BGCVectorDouble thermal_conductivity;
    BGCVectorDouble volume;
    BGCVectorDouble depth;
    BGCVectorDouble dz;
    BGCVectorDouble plant_wilting_factor;
    BGCVectorDouble rooting_depth_fraction;
    double shortwave_radiation;
    double longwave_radiation;
    double air_temperature;
    double vapor_pressure_air;
    double wind_speed;
    double precipitation;
    double elevation;
    double aspect;
    double slope;
    double atm_n2;
    double atm_o2;
    double atm_co2;
    double atm_ch4;
    double atm_n2o;
    double atm_h2;
    double atm_nh3;
  } BGCProperties;
  */

  typedef struct {
    BGCVectorInt aux_ints;  /* [-] */
    BGCVectorDouble aux_doubles;  /* [-] */
  } BGCAuxiliaryData;

  typedef struct {
    /* read data files/structures, initialize memory, basis management
       (includes reading database, swapping basis, etc.) */
    void (*DataTest)();

    void (*Setup)(
      BGCProperties* properties,
      BGCState* state,
      BGCSizes* sizes,
      int num_iterations,
      int num_columns,
      int ncells_per_col_);

    /* gracefully shutdown the engine, cleanup memory */
    void (*Shutdown)();

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
      double delta_t,
      BGCProperties* properties,
      BGCState* state,
      BGCSizes* sizes,
      int num_iterations,
      int num_columns);

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

  void CreateBGCInterface(const char* const engine_name, BGCInterface* interface);
  /*
  typedef struct {
    int error;
    char* message;
    bool converged;
    int num_rhs_evaluations;
    int num_jacobian_evaluations;
    int num_newton_iterations;
  } BGCEngineStatus;

  typedef struct {
    bool thread_safe;
    bool temperature_dependent;
    bool pressure_dependent;
    bool porosity_update;
    bool operator_splitting;
    bool global_implicit;
    int index_base;
  } BGCEngineFunctionality;

  typedef struct {
    AlquimiaVectorString primary_names;
    AlquimiaVectorInt    positivity;
    AlquimiaVectorString mineral_names;
    AlquimiaVectorString surface_site_names;
    AlquimiaVectorString ion_exchange_names;
    AlquimiaVectorString isotherm_species_names;
    AlquimiaVectorString aqueous_kinetic_names;
  } BGCProblemMetaData;

  typedef struct {
    double pH;
    AlquimiaVectorDouble aqueous_kinetic_rate;  // [?]
    AlquimiaVectorDouble mineral_saturation_index;  // [mol/sec/m^3]
    AlquimiaVectorDouble mineral_reaction_rate;  // [mol/sec/m^3 bulk]
    AlquimiaVectorDouble primary_free_ion_concentration; // [molality]
    AlquimiaVectorDouble primary_activity_coeff; // [-]
    AlquimiaVectorDouble secondary_free_ion_concentration; // [molality]
    AlquimiaVectorDouble secondary_activity_coeff; // [-]
  } BGCAuxiliaryOutputData;
  */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* ALQUIMIA_CONTAINERS_H_ */
