# -*- mode: cmake -*-
#
#  ATS
#    Constitutive relations for EcoSIM
#

# collect all sources
set(ats_ecosim_relations_src_files
  bulk_density/bulk_density_evaluator.cc
  bulk_density/bulk_density_model.cc
  hydraulic_conductivity/hydraulic_conductivity_evaluator.cc
  hydraulic_conductivity/hydraulic_conductivity_model.cc
)

set(ats_ecosim_relations_inc_files
  bulk_density/bulk_density_evaluator.hh
  bulk_density/bulk_density_model.hh
  hydraulic_conductivity/hydraulic_conductivity_evaluator.hh
  hydraulic_conductivity/hydraulic_conductivity_model.hh
)

set(ats_ecosim_relations_link_libs
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
add_amanzi_library(ats_ecosim_relations
                   SOURCE ${ats_ecosim_relations_src_files}
                   HEADERS ${ats_ecosim_relations_inc_files}
		   LINK_LIBS ${ats_ecosim_relations_link_libs})

register_evaluator_with_factory(
  HEADERFILE constitutive_models/bulk_density/bulk_density_evaluator_reg.hh
  LISTNAME ATS_ECOSIM_RELATIONS_REG
)

register_evaluator_with_factory(
  HEADERFILE constitutive_models/hydraulic_conductivity/hydraulic_conductivity_evaluator_reg.hh
  LISTNAME ATS_ECOSIM_RELATIONS_REG
)

generate_evaluators_registration_header(
    HEADERFILE ats_ecosim_relations_registration.hh
    LISTNAME   ATS_ECOSIM_RELATIONS_REG
    INSTALL    True
    )
