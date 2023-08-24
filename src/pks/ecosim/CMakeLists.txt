# -*- mode: cmake -*-

add_subdirectory(constitutive_relations)
add_subdirectory(data)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/data)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/bulk_density)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim/constitutive_relations/hydraulic_conductivity)
#include_directories(${TPL_INSTALL_PREFIX}/ecosim/f90src/ATSUtils/)

include_directories(${ECOSIM_INCLUDE_DIRS})

message(eco_inc_file="${ECOSIM_INCLUDE_DIRS}")

# ATS EcoSIM pk
# For adding F90 files add the F90 file to the ecosim source files and
# inc files. The inc file also needs a header

set(ats_ecosim_src_files
    EcoSIM_ATS_interface.cc
    BGCEngine.cc
    ecosim_wrappers.F90
    data/bgc_fortran_memory_mod.F90
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
