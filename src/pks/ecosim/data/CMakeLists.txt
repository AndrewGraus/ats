# -*- mode: cmake -*-
#
#  ATS
#    Data management programs for EcoSIM_ATS
#

# collect all sources
#set(ats_ecosim_data_src_files
#    BGC_constants.cc
#    BGC_memory.cc
#    BGC_containers.cc
#    bgc_fortran_memory_mod.F90
#)

set(TPL_INSTALL_PREFIX "/global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim")

include_directories(${TPL_INSTALL_PREFIX}/f90src/ATSUtils/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/utils/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/minimath/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelconfig/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/mesh/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelpars/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_datatype/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SoilPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SurPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_mods/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/PhysData/)

file(GLOB ats_ecosim_data_src_files
    BGC_constants.cc
    BGC_memory.cc
    BGC_containers.cc
    bgc_fortran_memory_mod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/*
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/minimath/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelconfig/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/mesh/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelpars/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/HydroTherm/SoilPhys/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/HydroTherm/SurPhys/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/HydroTherm/PhysData/*.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_mods/*.F90
)

set(ats_ecosim_data_inc_files
  BGC_constants.hh
  BGC_containers.hh
  BGC_memory.hh
)

set(ats_ecosim_data_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  operators
  solvers
  time_integration
  state
  pks
  chemistry_pk
  ats_pks
  ats_eos
  ats_operators
)

# make the library
add_amanzi_library(ats_ecosim_data
                   SOURCE ${ats_ecosim_data_src_files}
                   HEADERS ${ats_ecosim_data_inc_files}
		   LINK_LIBS ${ats_ecosim_data_link_libs})

generate_evaluators_registration_header(
    HEADERFILE ats_ecosim_data_registration.hh
    LISTNAME   ATS_ECOSIM_DATA_REG
    INSTALL    True
    )
