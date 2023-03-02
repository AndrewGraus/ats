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
    int num_primary;
    int num_sorbed;
    int num_minerals;
    int num_aqueous_complexes;
    int num_aqueous_kinetics;
    int num_surface_sites;
    int num_ion_exchange_sites;
    int num_isotherm_species;
    int num_aux_integers;
    int num_aux_doubles;
  } BGCSizes;

  typedef struct {
    BGCVectorDouble fluid_density;
    BGCVectorDouble gas_density;
    BGCVectorDouble ice_density;
    BGCVectorDouble porosity;
    BGCVectorDouble water_content;
    BGCVectorDouble temperature;
    BGCVectorDouble total_mobile;
  } BGCState;

  typedef struct {
    BGCVectorDouble liquid_saturation;
    BGCVectorDouble gas_saturation;
    BGCVectorDouble ice_saturation;
    BGCVectorDouble elevation;
    BGCVectorDouble relative_permeability;
    BGCVectorDouble conductivity;
    BGCVectorDouble volume;
  } BGCProperties;

  typedef struct {
    AlquimiaVectorInt aux_ints;  /* [-] */
    AlquimiaVectorDouble aux_doubles;  /* [-] */
  } BGCAuxiliaryData;
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
