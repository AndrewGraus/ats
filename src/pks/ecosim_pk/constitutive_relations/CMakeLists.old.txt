# -*- mode: cmake -*-
#
#  ATS
#    Constitutive relations for flow
#


# collect all sources

list(APPEND subdirs bulk_density hydraulic_conductivity)

set(ats_ecosim_relations_src_files "")
set(ats_ecosim_relations_inc_files "")

foreach(lcv IN LISTS subdirs)
  include_directories(${ATS_SOURCE_DIR}/src/pks/ecosim_pk/constitutive_relations/${lcv})

  file(GLOB subdir_sources "./${lcv}/*.cc")
  set(ats_ecosim_relations_src_files ${ats_ecosim_relations_src_files} ${subdir_sources})

  file(GLOB subdir_incs "./${lcv}/*.hh")
  set(ats_ecosim_relations_inc_files ${ats_ecosim_relations_inc_files} ${subdir_incs})

  file(GLOB registrations "./${lcv}/*_reg.hh" )
  foreach(reg_lcv IN LISTS registrations)
    register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME ATS_ECOISIM_RELATIONS_REG)
  endforeach(reg_lcv)

endforeach(lcv)

set(ats_ecosim_relations_link_libs
  ${Teuchos_LIBRARIES}
  ${Epetra_LIBRARIES}
  error_handling
  atk
  mesh
  data_structures
  whetstone
  solvers
  state
  )

# make the library
add_amanzi_library(ats_ecosim_relations
                   SOURCE ${ats_ecosim_relations_src_files}
                   HEADERS ${ats_ecosim_relations_inc_files}
		   LINK_LIBS ${ats_ecosim_relations_link_libs})

generate_evaluators_registration_header(
    HEADERFILE ats_ecosim_relations_registration.hh
    LISTNAME   ATS_ECOSIM_RELATIONS_REG
    INSTALL    True
    )
