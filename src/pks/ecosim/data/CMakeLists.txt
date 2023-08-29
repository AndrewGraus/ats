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
set(ECOSIM_LIB_LOCATION "/global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/local/lib")

include_directories(${TPL_INSTALL_PREFIX}/f90src/ATSUtils/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/utils/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/minimath/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelconfig/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/mesh/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelpars/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_datatype/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SoilPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SurfPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_mods/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/PhysData/)

set(ats_ecosim_data_src_files
    BGC_constants.cc
    BGC_memory.cc
    BGC_containers.cc
    bgc_fortran_memory_mod.F90
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
  ${ECOSIM_LIB_LOCATION}/libAPIData.a
  ${ECOSIM_LIB_LOCATION}/libAPIs.a
  ${ECOSIM_LIB_LOCATION}/libATSUtils_mods.a
  ${ECOSIM_LIB_LOCATION}/libPhysData.a
  ${ECOSIM_LIB_LOCATION}/libSnowPhys.a
  ${ECOSIM_LIB_LOCATION}/libSoilPhys.a
  ${ECOSIM_LIB_LOCATION}/libSurfPhys.a
  ${ECOSIM_LIB_LOCATION}/libbalances.a
  ${ECOSIM_LIB_LOCATION}/libbox_chem.a
  ${ECOSIM_LIB_LOCATION}/libbox_model.a
  ${ECOSIM_LIB_LOCATION}/libboxshared.a
  ${ECOSIM_LIB_LOCATION}/libdisturb_mods.a
  ${ECOSIM_LIB_LOCATION}/libecosim_datatype.a
  ${ECOSIM_LIB_LOCATION}/libecosim_main.a
  ${ECOSIM_LIB_LOCATION}/libecosim_mods.a
  ${ECOSIM_LIB_LOCATION}/libioutils.a
  ${ECOSIM_LIB_LOCATION}/liblayers_chem.a
  ${ECOSIM_LIB_LOCATION}/liblayers_model.a
  ${ECOSIM_LIB_LOCATION}/libmesh.a
  ${ECOSIM_LIB_LOCATION}/libminimath.a
  ${ECOSIM_LIB_LOCATION}/libmodelconfig.a
  ${ECOSIM_LIB_LOCATION}/libmodelforc.a
  ${ECOSIM_LIB_LOCATION}/libmodelpars.a
  ${ECOSIM_LIB_LOCATION}/libnonsalt.a
  ${ECOSIM_LIB_LOCATION}/libplant_bgc.a
  ${ECOSIM_LIB_LOCATION}/libplant_driver.a
  ${ECOSIM_LIB_LOCATION}/libsalt.a
  ${ECOSIM_LIB_LOCATION}/libutil.a
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
