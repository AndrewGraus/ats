# -*- mode: cmake -*-

add_subdirectory(constitutive_relations)

include_directories(${ATS_SOURCE_DIR}/src/pks)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/root_storage)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/flux_uptake)
include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/atm_flux)

# ATS EcoSIM pk

set(ats_ecosim_src_files
    Ecosim_interface_test.cc
    constitutive_relations/flux_uptake/eco_fort_test.F90
  )

set(ats_ecosim_inc_files
    Ecosim_interface_test.hh
    constitutive_relations/flux_uptake/eco_fort_test.F90
    constitutive_relations/flux_uptake/eco_fort_test_interface_private.hh
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
  ecosim_relations
  )


add_amanzi_library(ats_ecosim
                   SOURCE ${ats_ecosim_src_files}
                   HEADERS ${ats_ecosim_inc_files}
		   LINK_LIBS ${ats_ecosim_link_libs})

#================================================
# register evaluators/factories/pks

register_evaluator_with_factory(
  HEADERFILE Ecosim_interface_test_reg.hh
  LISTNAME ATS_ECOSIM_REG
)

#Do I need this?
#seems like you don't need the file, it's auto-generated
generate_evaluators_registration_header(
  HEADERFILE ats_ecosim_registration.hh
  LISTNAME   ATS_ECOSIM_REG
  INSTALL    True
)
