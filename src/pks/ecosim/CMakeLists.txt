# -*- mode: cmake -*-

add_subdirectory(constitutive_relations)
add_subdirectory(data)

set(ECOSIM_INSTALL_PREFIX "/global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim")
set(ECOSIM_LIB_LOCATION "/global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/local/lib")

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/data)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/bulk_density)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/hydraulic_conductivity)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ATSUtils/)

include_directories(${TPL_INSTALL_PREFIX}/f90src/utils/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/minimath/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelconfig/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/mesh/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/modelpars/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_datatype/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SoilPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/SurfPhys/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/HydroTherm/PhysData/)
include_directories(${TPL_INSTALL_PREFIX}/f90src/ecosim_mods/)

include_directories(${ECOSIM_INCLUDE_DIRS})

message(eco_inc_file="${ECOSIM_INCLUDE_DIRS}")

# ATS EcoSIM pk
# For adding F90 files add the F90 file to the ecosim source files and
# inc files. The inc file also needs a header

#set(ats_ecosim_src_files
#    EcoSIM_ATS_interface.cc
#    BGCEngine.cc
#    ecosim_wrappers.F90
#    data/bgc_fortran_memory_mod.F90
#)

set(ats_ecosim_src_files
    EcoSIM_ATS_interface.cc
    BGCEngine.cc
    ecosim_wrappers.F90
    data/bgc_fortran_memory_mod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/BGC_containers.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSCPLMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/ATSEcoSIMInitMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/SharedDataMod.F90
    ${ECOSIM_INSTALL_PREFIX}/f90src/ATSUtils/c_f_interface_module.F90
)

set(ats_ecosim_inc_files
    ecosim_mod_test_wrapper.h
    EcoSIM_ATS_interface.hh
    BGCEngine.hh
    ecosim_interface.h
  )

set(ats_ecosim_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  ${ECOSIM_LIBRARIES}
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
  ats_ecosim_data
  ats_ecosim_relations
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

message(inc_files="${ats_ecosim_inc_files}")

add_amanzi_library(ats_ecosim
                   SOURCE ${ats_ecosim_src_files}
                   HEADERS ${ats_ecosim_inc_files}
		   LINK_LIBS ${ats_ecosim_link_libs})

#================================================
# register evaluators/factories/pks

register_evaluator_with_factory(
  HEADERFILE EcoSIM_ATS_interface_reg.hh
  LISTNAME ATS_ECOSIM_REG
)

generate_evaluators_registration_header(
  HEADERFILE ats_ecosim_registration.hh
  LISTNAME   ATS_ECOSIM_REG
  INSTALL    True
)
