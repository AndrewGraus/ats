# -*- mode: cmake -*-
#
#  ATS
#    Data management programs for EcoSIM_ATS
#

# collect all sources
set(ats_ecosim_data_src_files
  BGC_constants.cc
  BGC_memory.cc
  BGC_containers.cc
  bgc_fortran_memory_mod.F90
  /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/BGC_containers.F90
  /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/ATSCPLMod.F90
  /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/c_f_interface_module.F90
  /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/data_kind_mod.F90
  /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/SharedDataMod.F90
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
