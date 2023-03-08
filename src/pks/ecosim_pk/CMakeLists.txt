# -*- mode: cmake -*-

add_subdirectory(constitutive_relations)
add_subdirectory(data)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/data)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/bulk_density)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/hydraulic_conductivity)

# ATS EcoSIM pk
# For adding F90 files add the F90 file to the ecosim source files and
# inc files. The inc file also needs a header

set(ats_ecosim_src_files
    EcoSIM_ATS_interface.cc
    BGCEngine.cc
  )

set(ats_ecosim_inc_files
    EcoSIM_ATS_interface.hh
    BGCEngine.hh
  )

set(ats_ecosim_link_libs
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
  ats_ecosim_data
  ats_ecosim_relations
  )


add_amanzi_library(ats_ecosim
                   SOURCE ${ats_ecosim_src_files}
                   HEADERS ${ats_ecosim_inc_files}
		   LINK_LIBS ${ats_ecosim_link_libs})

#================================================
# register evaluators/factories/pks

register_evaluator_with_factory(
  HEADERFILE ats_ecosim_data_registration.hh
  LISTNAME ATS_ECOSIM_DATA_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/bulk_density/bulk_density_evaluator_reg.hh
  LISTNAME ATS_ECOSIM_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_relations/hydraulic_conductivity/hydraulic_conductivity_evaluator_reg.hh
  LISTNAME ATS_ECOSIM_REG
)

register_evaluator_with_factory(
  HEADERFILE EcoSIM_interface_reg.hh
  LISTNAME ATS_ECOSIM_REG
)

generate_evaluators_registration_header(
  HEADERFILE ats_ecosim_registration.hh
  LISTNAME   ATS_ECOSIM_REG
  INSTALL    True
)
