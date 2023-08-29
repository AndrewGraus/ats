# -*- mode: cmake -*-

add_subdirectory(constitutive_relations)
add_subdirectory(data)

set(TPL_INSTALL_PREFIX "/global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim")

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
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/ATSCPLMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/ATSEcoSIMInitMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/BGC_containers.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ATSUtils/SharedDataMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/data_kind_mod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/abortutils.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/fileUtil.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/ecosim_time_mod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/EcoSimConst.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/ncdio_pio.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/utils/ecosim_log_mod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/SoilWaterDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/GridDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/SOMDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/SoilPhysDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/LandSurfDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/ClimForcDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_datatype/SoilPropertyDataType.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/ecosim_mods/InitAllocMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/mesh/GridConsts.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelconfig/EcoSIMConfig.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelconfig/TracerIDMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelconfig/ElmIDMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelconfig/EcoSIMCtrlMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/minimath/MiniMathMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/GrosubsMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/PhotoSynsMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/RootMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/LitterFallMod.f90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/PlantBranchMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelpars/EcoSiMParDataMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/modelpars/GrosubPars.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/plant_bgc/NoduleBGCMod.F90
    /global/home/users/agraus/code/ats_dev_dir/amanzi_tpls-install-master-Debug/ecosim/f90src/APIs/PlantAPIData.F90
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
