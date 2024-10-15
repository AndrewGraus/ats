/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus

  This is the main PK for the EcoSIM-ATS interface. This code is written
  following the example of Alquimia with some additional code from the
  SimpleBGC code for walking the columns.

  The idea is to take the basic code used by alquimia and repurpose it so
  that it works on a column by column basis instead of a cell by cell basis

  --------------------------------------------------------------------------*/

#include <algorithm>
#include <set>
#include <string>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "Mesh.hh"

// include custom evaluators here
#include "hydraulic_conductivity_evaluator.hh"
#include "bulk_density_evaluator.hh"

#include "pk_helpers.hh"
#include "EcoSIM_ATS_interface.hh"

namespace Amanzi {
namespace EcoSIM {

EcoSIM::EcoSIM(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution):
  PK_Physical(pk_tree, global_list, S, solution),
  PK(pk_tree, global_list, S, solution),
  ncells_per_col_(-1),
  saved_time_(0.0)
  {
    //grab the surface and subsurface domains
    domain_ = plist_->get<std::string>("domain name", "domain");
    domain_surface_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

    // transport
    tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");
    //Remember tcc components are accessed by tcc[i][c] where i is the component and c is the cell

    //Flow
    porosity_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_liquid_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    saturation_gas_key_ = Keys::readKey(*plist_,domain_,"saturation gas", "saturation_gas");
    saturation_ice_key_ = Keys::readKey(*plist_,domain_,"saturation ice", "saturation_ice");
    water_content_key_ = Keys::readKey(*plist_,domain_,"water content","water_content");
    relative_permeability_key_ = Keys::readKey(*plist_,domain_,"relative permeability","relative_permeability");
    matric_pressure_key_ = Keys::readKey(*plist_,domain_,"matric pressure","matric_pressure");
    liquid_density_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
    ice_density_key_ = Keys::readKey(*plist_, domain_, "mass density ice", "mass_density_ice");
    gas_density_key_ = Keys::readKey(*plist_, domain_,"mass density gas", "mass_density_gas");
    gas_density_key_test_ = Keys::readKey(*plist_, domain_, "mass density gas", "mass_density_gas");
    rock_density_key_ = Keys::readKey(*plist_, domain_, "density rock", "density_rock");

    //energy
    T_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
    thermal_conductivity_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");

    //Sources
    surface_water_source_key_ = Keys::readKey(*plist_, domain_surface_, "surface water source", "water_source");
    surface_energy_source_key_ =
      Keys::readKey(*plist_, domain_surface_, "surface energy source", "total_energy_source");
    subsurface_water_source_key_ =
      Keys::readKey(*plist_, domain_, "subsurface water source", "water_source");
    subsurface_energy_source_key_ =
      Keys::readKey(*plist_, domain_, "subsurface energy source", "total_energy_source");
    surface_energy_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_surface_, "surface energy source ecosim", "ecosim_source");
    surface_water_source_ecosim_key_ =
      Keys::readKey(*plist_, domain_surface_, "surface water source ecosim", "ecosim_water_source");

    //Snow

    //Other
    cell_volume_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
    //ecosim_aux_data_key_ = Keys::readKey(*plist_, domain_, "ecosim aux data", "ecosim_aux_data");
    f_wp_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    f_root_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");

    //Custom Evaluator keys
    hydraulic_conductivity_key_ = Keys::readKey(*plist_, domain_, "hydraulic conductivity", "hydraulic_conductivity");
    bulk_density_key_ = Keys::readKey(*plist_, domain_, "bulk density", "bulk_density");

    //Surface balance items
    sw_key_ =
      Keys::readKey(*plist_, domain_surface_, "incoming shortwave radiation", "incoming_shortwave_radiation");
    lw_key_ =
      Keys::readKey(*plist_,domain_surface_, "incoming longwave radiation", "incoming_longwave_radiation");
    air_temp_key_ = Keys::readKey(*plist_, domain_surface_, "air temperature", "air_temperature");
    vp_air_key_ = Keys::readKey(*plist_, domain_surface_, "vapor pressure air", "vapor_pressure_air");
    wind_speed_key_ = Keys::readKey(*plist_, domain_surface_, "wind speed", "wind_speed");
    p_rain_key_ = Keys::readKey(*plist_, domain_surface_, "precipitation rain", "precipitation_rain");
    elev_key_ = Keys::readKey(*plist_, domain_surface_, "elevation", "elevation");
    aspect_key_ = Keys::readKey(*plist_, domain_surface_, "aspect", "aspect");
    slope_key_ = Keys::readKey(*plist_, domain_surface_, "slope", "slope_magnitude");
    snow_depth_key_ = Keys::readKey(*plist_, domain_surface_, "snow depth", "snow_depth");

    //Atmospheric abundance keys
    atm_n2_ = plist_->get<double>("atmospheric N2");
    atm_o2_ = plist_->get<double>("atmospheric O2");
    atm_co2_ = plist_->get<double>("atmospheric CO2");
    atm_ch4_ = plist_->get<double>("atmospheric CH4");
    atm_n2o_ = plist_->get<double>("atmospheric N2O");
    atm_h2_ = plist_->get<double>("atmospheric H2");
    atm_nh3_ = plist_->get<double>("atmospheric NH3");
    pressure_at_field_capacity = plist_->get<double>("Field Capacity [Mpa]");
    pressure_at_wilting_point = plist_->get<double>("Wilting Point [Mpa]");

    dt_ = plist_->get<double>("initial time step");
    c_m_ = plist_->get<double>("heat capacity [MJ mol^-1 K^-1]");

    //This initialized the engine (found in BGCEngine.cc) This is the code that
    //actually points to the driver
    if (!plist_->isParameter("engine")) {
      Errors::Message msg;
      msg << "No 'engine' parameter found in the parameter list for 'BGC'.\n";
      Exceptions::amanzi_throw(msg);
    }
    if (!plist_->isParameter("engine input file")) {
      Errors::Message msg;
      msg << "No 'engine input file' parameter found in the parameter list for 'BGC'.\n";
      Exceptions::amanzi_throw(msg);
    }
    std::string engine_name = plist_->get<std::string>("engine");
    std::string engine_inputfile = plist_->get<std::string>("engine input file");
    bgc_engine_ = Teuchos::rcp(new BGCEngine(engine_name, engine_inputfile));

    /*(AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
    Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
    wrms_ = createWRMPartition(sublist);
    InitializeFromPlist_();
    */
  }


// -- Destroy ansilary data structures.
EcoSIM::~EcoSIM()
  {
  if (bgc_initialized_)
    bgc_engine_->FreeState(bgc_props_, bgc_state_, bgc_aux_data_);
  }

// -- Setup step
void EcoSIM::Setup() {
  //Need to do some basic setup of the columns:
  mesh_surf_ = S_->GetMesh(domain_surface_);
  num_columns_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (unsigned int column = 0; column != num_columns_; ++column) {
    int f = mesh_surf_->entity_get_parent(AmanziMesh::CELL, column);
    auto& col_iter = mesh_->cells_of_column(column);
    std::size_t ncol_cells = col_iter.size();

    double column_area = mesh_->face_area(f);

    if (ncells_per_col_ < 0) {
      ncells_per_col_ = ncol_cells;
    } else {
      AMANZI_ASSERT(ncol_cells == ncells_per_col_);
    }
  }

  //This is for the Auxiliary Data which we will need
  //commenting for now because we don't need it yet
  /*
  if (plist_->isParameter("auxiliary data")) {
    auto names = plist_->get<Teuchos::Array<std::string> >("auxiliary data");

    for (auto it = names.begin(); it != names.end(); ++it) {
      Key aux_field_name = Keys::getKey(domain_, *it);
      if (!S_->HasRecord(aux_field_name)) {
        S_->Require<CompositeVector, CompositeVectorSpace>(aux_field_name, tag_next_, passwd_)
          .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, 1);
      }
    }
  }

  // Setup more auxiliary data
  if (!S_->HasRecord(alquimia_aux_data_key_, tag_next_)) {
    int num_aux_data = bgc_engine_->Sizes().num_aux_integers + bgc_engine_->Sizes().num_aux_doubles;
    S_->Require<CompositeVector, CompositeVectorSpace>(alquimia_aux_data_key_, tag_next_, passwd_)
      .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, num_aux_data);

    S_->GetRecordW(alquimia_aux_data_key_, tag_next_, passwd_).set_io_vis(false);
  }
  */

  //Setup Evaluators
  requireAtNext(hydraulic_conductivity_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtCurrent(hydraulic_conductivity_key_, tag_current_, *S_, name_);

  requireAtNext(bulk_density_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtCurrent(bulk_density_key_, tag_current_, *S_, name_);

  requireAtNext(matric_pressure_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtCurrent(matric_pressure_key_, tag_current_, *S_, name_);

  /*S_->RequireEvaluator(test_key_, tag_current_);
  S_->Require<CompositeVector, CompositeVectorSpace>(test_key_, tag_current_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::CELL, 1);  
  S_->RequireEvaluator(test_key_, tag_current_);
  */

  S_->RequireEvaluator(snow_depth_key_, tag_current_);
  S_->Require<CompositeVector, CompositeVectorSpace>(snow_depth_key_, tag_current_)
    .SetMesh(mesh_surf_)
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(snow_depth_key_, tag_current_);  

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Setup of PK was successful"
    << vo_->reset() << std::endl << std::endl;
  }
}

// -- Initialize owned (dependent) variables.
void EcoSIM::Initialize() {

  //Need to know the number of components to initialize data structures
  const Epetra_MultiVector& tcc= *(S_->GetPtr<CompositeVector>(tcc_key_, Tags::DEFAULT)->ViewComponent("cell"));
  int tcc_num = tcc.NumVectors();
  Teuchos::OSTab tab = vo_->getOSTab();

  num_columns_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //Now we call the engine's init state function which allocates the data
  bgc_engine_->InitState(bgc_props_, bgc_state_, bgc_aux_data_, ncells_per_col_, tcc_num, num_columns_);

  int ierr = 0;

  /*if (S_->HasRecord(suc_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "has suction key." << std::endl;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Does not have suction key at default." << std::endl;
  }
  if (S_->HasRecord(suc_key_, Tags::CURRENT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "has suction key at current" << std::endl;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Does not have suction key at current" << std::endl;
  }
  if (S_->HasRecord(suc_key_, Tags::NEXT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "has suction key at next" << std::endl;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Does not have suction key at next" << std::endl;
  }*/

  // Ensure dependencies are filled
  // May not need to update (also causes an assertion error if called before
  // The PK that owns the variable
  /*S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(porosity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(relative_permeability_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_wp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_root_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(subsurface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(subsurface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  */
  //S_->GetEvaluator(suc_key_, Tags::DEFAULT).Update(*S_, name_);

  //Surface properties from met data
  /*S_->GetEvaluator(sw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(lw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(p_rain_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(air_temp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(vp_air_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(wind_speed_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(aspect_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(slope_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(surface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(surface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);

  if (S_->HasRecord(gas_density_key_test_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found mass density gas key" << std::endl;
    S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(gas_density_key_, Tags::DEFAULT).Update(*S_, name_);
    has_gas = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Did not find gas key" << std::endl;
    has_gas = false;
  }

  if (S_->HasRecord(T_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found temp key" << std::endl;
    S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(thermal_conductivity_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(cell_volume_key_, Tags::DEFAULT).Update(*S_, name_);
    has_energy = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Did not find temp key" << std::endl;
    has_energy = false;
  }
  */
  /*if (S_->HasRecord(surface_snow_depth_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found surface_snow_depth_key" << std::endl;
  } else if (S_->HasRecord(snow_depth_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found snow_depth_key" << std::endl;
  } else {
    *vo_->os() << "neither snow depth key found" << std::endl;
  }*/

  if (S_->HasRecord(ice_density_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found ice key" << std::endl;
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_density_key_, Tags::DEFAULT).Update(*S_, name_);
    has_ice = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Did not find ice key" << std::endl;
    has_ice = false;
  }

  //Initialize owned evaluators
  S_->GetW<CompositeVector>(hydraulic_conductivity_key_, Tags::DEFAULT, "hydraulic_conductivity").PutScalar(1.0);
  S_->GetRecordW(hydraulic_conductivity_key_, Tags::DEFAULT, "hydraulic_conductivity").set_initialized();

  S_->GetW<CompositeVector>(bulk_density_key_, Tags::DEFAULT, "bulk_density").PutScalar(1.0);
  S_->GetRecordW(bulk_density_key_, Tags::DEFAULT, "bulk_density").set_initialized();

  S_->GetW<CompositeVector>(matric_pressure_key_, Tags::DEFAULT, "matric_pressure").PutScalar(1.0);
  S_->GetRecordW(matric_pressure_key_, Tags::DEFAULT, "matric_pressure").set_initialized();

  int num_columns_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //Looping over the columns and initializing
  /*for (int col=0; col!=num_columns_; ++col) {
    ierr = InitializeSingleColumn(col);
  }*/

  //loop over processes instead:
  num_columns_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  num_columns_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  //Trying to loop over processors now:
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

      InitializeSingleProcess(p_rank);
    }
  }

  // verbose message
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T="
        << S_->get_time() << vo_->reset() << std::endl << std::endl;
  }
}

void EcoSIM::CommitStep(double t_old, double t_new, const Tag& tag) {

  // I don't know that we will have much to do here. In SimpleBGC they just copy
  // Data to the pfts, which we won't be doing. In Alquimia they just save the time
  // As below.

  saved_time_ = t_new;

}

bool EcoSIM::AdvanceStep(double t_old, double t_new, bool reinit) {
  double dt = t_new - t_old;
  current_time_ = saved_time_ + dt;

  // If we are given a dt that is less than the one we wanted, we don't record it.
  // This is in alquimia but I'm not quite sure why
  //if (dt < dt_next_) {
  //  dt_prev_ = dt_next_;
  //} else {
  //  dt_prev_ = dt;
  //}

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  
  // Ensure dependencies are filled
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(porosity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(relative_permeability_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_density_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cell_volume_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_wp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_root_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(subsurface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(subsurface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(matric_pressure_key_, Tags::DEFAULT).Update(*S_, name_);

  //Surface data from met data
  S_->GetEvaluator(sw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(lw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(air_temp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(vp_air_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(wind_speed_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(p_rain_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(aspect_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(slope_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(surface_energy_source_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(surface_water_source_key_, Tags::DEFAULT).Update(*S_, name_);

  if (has_gas) {
    S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(gas_density_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  if (has_ice) {
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_density_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(thermal_conductivity_key_, Tags::DEFAULT).Update(*S_, name_);

  //Update owned evaluators
  Teuchos::RCP<const CompositeVector> hydra_cond = S_->GetPtr<CompositeVector>(hydraulic_conductivity_key_, Tags::DEFAULT);
  S_->GetEvaluator(hydraulic_conductivity_key_, Tags::DEFAULT).Update(*S_, name_);

  Teuchos::RCP<const CompositeVector> bulk_dens = S_->GetPtr<CompositeVector>(bulk_density_key_, Tags::DEFAULT);
  S_->GetEvaluator(bulk_density_key_, Tags::DEFAULT).Update(*S_, name_);

  Teuchos::RCP<const CompositeVector> matric_pressure = S_->GetPtr<CompositeVector>(matric_pressure_key_, Tags::DEFAULT);
  S_->GetEvaluator(matric_pressure_key_, Tags::DEFAULT).Update(*S_, name_);

  AmanziMesh::Entity_ID num_columns_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // grab the required fields
  S_->GetEvaluator("porosity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& porosity = *(*S_->Get<CompositeVector>("porosity", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_saturation = *(*S_->Get<CompositeVector>("saturation_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("water_content", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& water_content = *(*S_->Get<CompositeVector>("water_content", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("relative_permeability", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& relative_permeability = *(*S_->Get<CompositeVector>("relative_permeability", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("matric_pressure", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& suction_head = *(*S_->Get<CompositeVector>("matric_pressure", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("mass_density_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_density = *(*S_->Get<CompositeVector>("mass_density_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("density_rock", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& rock_density = *(*S_->Get<CompositeVector>("density_rock", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& cell_volume = *(*S_->Get<CompositeVector>("cell_volume", tag_next_)
          .ViewComponent("cell",false))(0);

  /*S_->GetEvaluator("plant_wilting_factor", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& plant_wilting_factor = *(*S_->Get<CompositeVector>("plant_wilting_factor", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("rooting_depth_fraction", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& rooting_depth_fraction = *(*S_->Get<CompositeVector>("rooting_depth_fraction", tag_next_)
          .ViewComponent("cell",false))(0);*/

  if (has_gas) {
    S_->GetEvaluator("mass_density_gas", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& gas_density = *(*S_->Get<CompositeVector>("mass_density_gas", tag_next_)
            .ViewComponent("cell",false))(0);

    S_->GetEvaluator("saturation_gas", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& gas_saturation = *(*S_->Get<CompositeVector>("saturation_gas", tag_next_)
            .ViewComponent("cell",false))(0);
  }

  //Atm abundances
  S_->GetEvaluator("surface-incoming_shortwave_radiation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& sw_rad = *(*S_->Get<CompositeVector>("surface-incoming_shortwave_radiation", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-incoming_longwave_radiation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& lw_rad = *(*S_->Get<CompositeVector>("surface-incoming_longwave_radiation", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-air_temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& t_air = *(*S_->Get<CompositeVector>("surface-air_temperature", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-vapor_pressure_air", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& p_vap = *(*S_->Get<CompositeVector>("surface-vapor_pressure_air", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-wind_speed", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& v_wind = *(*S_->Get<CompositeVector>("surface-wind_speed", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-precipitation_rain", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& p_rain = *(*S_->Get<CompositeVector>("surface-precipitation_rain", tag_next_)
          .ViewComponent("cell",false))(0);

  //S_->GetEvaluator("surface-precipitation_snow", tag_next_).Update(*S_, name_);
  //const Epetra_MultiVector& p_snow = *(*S_->Get<CompositeVector>("surface-precipitation_snow", tag_next_)
  //        .ViewComponent("cell",false))(0);

  S_->GetEvaluator("surface-elevation", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& elevation = *S_->Get<CompositeVector>("surface-elevation", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-aspect", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& aspect = *S_->Get<CompositeVector>("surface-aspect", tag_next_)
          .ViewComponent("cell",false);

  S_->GetEvaluator("surface-slope_magnitude", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& slope = *S_->Get<CompositeVector>("surface-slope_magnitude", tag_next_)
          .ViewComponent("cell",false);

  if (has_ice) {
    S_->GetEvaluator("mass_density_ice", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& ice_density = *(*S_->Get<CompositeVector>("mass_density_ice", tag_next_)
            .ViewComponent("cell",false))(0);

    S_->GetEvaluator("saturation_ice", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& ice_saturation = *(*S_->Get<CompositeVector>("saturation_ice", tag_next_)
            .ViewComponent("cell",false))(0);
  }

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& temp = *(*S_->Get<CompositeVector>("temperature", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("thermal_conductivity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& thermal_conductivity = *(*S_->Get<CompositeVector>("thermal_conductivity", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator(snow_depth_key_, tag_current_).Update(*S_, name_);
  const Epetra_MultiVector& snow_depth =
  *S_->Get<CompositeVector>(snow_depth_key_, tag_current_).ViewComponent("cell", false);

  //loop over processes instead:
  num_columns_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  num_columns_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  //Trying to loop over processors now:
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

      AdvanceSingleProcess(dt, p_rank);
    }
  }

}

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_Vector& vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(column);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = vec[vec_index];
  }
}

void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID column, const Teuchos::Ptr<Epetra_SerialDenseVector> vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(column);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = (*vec)[vec_index];
  }
}

void EcoSIM::MatrixFieldToColumn_(AmanziMesh::Entity_ID column, const Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr)
  {
    int n_comp = m_arr.NumVectors();
    auto& col_iter = mesh_->cells_of_column(column);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        (*col_arr)(i,j) = m_arr[j][col_iter[i]];
      }
    }
  }

// helper function for pushing column back to field
void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID column, Epetra_Vector& vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(column);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    vec[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID column, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(column);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*vec)[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::MatrixColumnToField_(AmanziMesh::Entity_ID column, Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr) {

    int n_comp = m_arr.NumVectors();
    auto& col_iter = mesh_->cells_of_column(column);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        m_arr[j][col_iter[i]] = (*col_arr)(i,j);
      }
    }

  }

// helper function for collecting column dz and depth
void EcoSIM::ColDepthDz_(AmanziMesh::Entity_ID column,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, column);
  auto& col_iter = mesh_->cells_of_column(column);
  ncells_per_col_ = col_iter.size();

  AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  Teuchos::OSTab tab = vo_->getOSTab();

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

    // dz
    // -- find face_below
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    double vol = mesh_->cell_volume(col_iter[i]);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }
}

// helper function for collecting dz, depth, and volume for a given column
void EcoSIM::VolDepthDz_(AmanziMesh::Entity_ID column,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz,
			    Teuchos::Ptr<Epetra_SerialDenseVector> volume) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, column);
  auto& col_iter = mesh_->cells_of_column(column);
  ncells_per_col_ = col_iter.size();

  AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

    // dz
    // -- find face_below
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    //double vol = mesh_->cell_volume(col_iter[i]);
    (*volume)[i] = mesh_->cell_volume(col_iter[i]);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }
}

//Copy to EcoSIM
void EcoSIM::CopyToEcoSIM_process(int proc_rank,
                                 BGCProperties& props,
                                 BGCState& state,
                                 BGCAuxiliaryData& aux_data,
                               const Tag& water_tag)
{
  //This is the copy function for a loop over a single process instead of a single column
  //Fill state with ATS variables that are going to be changed by EcoSIM
  const Epetra_Vector& porosity = *(*S_->Get<CompositeVector>(porosity_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_MultiVector& tcc= *(S_->GetPtr<CompositeVector>(tcc_key_, water_tag)->ViewComponent("cell"));
  int tcc_num = tcc.NumVectors();

  const Epetra_Vector& liquid_saturation = *(*S_->Get<CompositeVector>(saturation_liquid_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& water_content = *(*S_->Get<CompositeVector>(water_content_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& relative_permeability = *(*S_->Get<CompositeVector>(relative_permeability_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_density = *(*S_->Get<CompositeVector>(liquid_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rock_density = *(*S_->Get<CompositeVector>(rock_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& cell_volume = *(*S_->Get<CompositeVector>(cell_volume_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& hydraulic_conductivity = *(*S_->Get<CompositeVector>(hydraulic_conductivity_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& matric_pressure = *(*S_->Get<CompositeVector>(matric_pressure_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& bulk_density = *(*S_->Get<CompositeVector>(bulk_density_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rooting_depth_fraction = *(*S_->Get<CompositeVector>(f_root_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& plant_wilting_factor = *(*S_->Get<CompositeVector>(f_wp_key_, water_tag).ViewComponent("cell", false))(0);

  const Epetra_Vector& shortwave_radiation = *(*S_->Get<CompositeVector>(sw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& longwave_radiation = *(*S_->Get<CompositeVector>(lw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& air_temperature = *(*S_->Get<CompositeVector>(air_temp_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& vapor_pressure_air = *(*S_->Get<CompositeVector>(vp_air_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& wind_speed= *(*S_->Get<CompositeVector>(wind_speed_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& precipitation = *(*S_->Get<CompositeVector>(p_rain_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& elevation = *(*S_->Get<CompositeVector>(elev_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& aspect = *(*S_->Get<CompositeVector>(aspect_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& slope = *(*S_->Get<CompositeVector>(slope_key_, water_tag).ViewComponent("cell", false))(0);

  const Epetra_Vector& surface_energy_source = *(*S_->Get<CompositeVector>(surface_energy_source_ecosim_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& subsurface_energy_source = *(*S_->Get<CompositeVector>(subsurface_energy_source_key_, water_tag).ViewComponent("cell", false))(0);

  const Epetra_Vector& surface_water_source = *(*S_->Get<CompositeVector>(surface_water_source_ecosim_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& subsurface_water_source = *(*S_->Get<CompositeVector>(subsurface_water_source_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& snow_depth = *(*S_->Get<CompositeVector>(snow_depth_key_, water_tag).ViewComponent("cell", false))(0);

  auto col_porosity = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_relative_permeability = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_mat_p = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_b_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_depth = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_dz = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rf = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_energy_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_water_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_depth_c = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));

  //Gather columns on this process:
  num_columns_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  num_columns_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  //Trying to loop over processors now:
  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //Loop over columns on this process
  for (int column=0; column!=num_columns_local; ++column) {
    FieldToColumn_(column,porosity,col_porosity.ptr());
    FieldToColumn_(column,liquid_saturation,col_l_sat.ptr());
    FieldToColumn_(column,water_content,col_wc.ptr());
    FieldToColumn_(column,relative_permeability,col_relative_permeability.ptr());
    FieldToColumn_(column,liquid_density,col_l_dens.ptr());
    FieldToColumn_(column,rock_density,col_r_dens.ptr());
    //FieldToColumn_(column,cell_volume,col_vol.ptr());
    FieldToColumn_(column,hydraulic_conductivity,col_h_cond.ptr());
    FieldToColumn_(column,bulk_density,col_b_dens.ptr());
    FieldToColumn_(column,plant_wilting_factor,col_wp.ptr());
    FieldToColumn_(column,rooting_depth_fraction,col_rf.ptr());
    FieldToColumn_(column,subsurface_water_source,col_ss_water_source.ptr());
    FieldToColumn_(column,subsurface_energy_source,col_ss_energy_source.ptr());
    FieldToColumn_(column,matric_pressure,col_mat_p.ptr());

    MatrixFieldToColumn_(column, tcc, col_tcc.ptr());

    if (has_gas) {
      const Epetra_Vector& gas_saturation = *(*S_->Get<CompositeVector>(saturation_gas_key_, water_tag).ViewComponent("cell", false))(0);
      const Epetra_Vector& gas_density = *(*S_->Get<CompositeVector>(gas_density_key_, water_tag).ViewComponent("cell", false))(0);

      FieldToColumn_(column,gas_saturation,col_g_sat.ptr());
      FieldToColumn_(column,gas_density,col_g_dens.ptr());
    }

    if (has_ice) {
      const Epetra_Vector& ice_saturation = *(*S_->Get<CompositeVector>(saturation_ice_key_, water_tag).ViewComponent("cell", false))(0);
      const Epetra_Vector& ice_density = *(*S_->Get<CompositeVector>(ice_density_key_, water_tag).ViewComponent("cell", false))(0);

      FieldToColumn_(column,ice_saturation,col_i_sat.ptr());
      FieldToColumn_(column,ice_density,col_i_dens.ptr());
    }

    const Epetra_Vector& temp = *(*S_->Get<CompositeVector>(T_key_, water_tag).ViewComponent("cell", false))(0);
    const Epetra_Vector& thermal_conductivity = *(*S_->Get<CompositeVector>(thermal_conductivity_key_, water_tag).ViewComponent("cell", false))(0);

    FieldToColumn_(column,temp, col_temp.ptr());
    FieldToColumn_(column,thermal_conductivity,col_cond.ptr());

    // This is for computing depth
    //ColDepthDz_(column, col_depth.ptr(), col_dz.ptr());

    VolDepthDz_(column, col_depth.ptr(), col_dz.ptr(), col_vol.ptr());
    double sum = 0.0;
    for (int i = ncells_per_col_ - 1; i >= 0; --i) {
        sum += (*col_dz)[i];
        (*col_depth_c)[i] = sum;
    }
   
    for (int i=0; i < ncells_per_col_; ++i) {
      state.liquid_density.data[column][i] = (*col_l_dens)[i];
      state.porosity.data[column][i] = (*col_porosity)[i];
      state.water_content.data[column][i] = (*col_wc)[i];
      state.hydraulic_conductivity.data[column][i] = (*col_h_cond)[i];
      state.bulk_density.data[column][i] = (*col_b_dens)[i];
      state.subsurface_water_source.data[column][i] = (*col_ss_water_source)[i];
      state.subsurface_energy_source.data[column][i] = (*col_ss_energy_source)[i];
      state.matric_pressure.data[column][i] = (*col_mat_p)[i];
      
      props.plant_wilting_factor.data[column][i] = (*col_wp)[i];
      props.rooting_depth_fraction.data[column][i] = (*col_rf)[i];
      props.liquid_saturation.data[column][i] = (*col_l_sat)[i];
      props.relative_permeability.data[column][i] = (*col_relative_permeability)[i];
      props.volume.data[column][i] = (*col_vol)[i];
      props.depth.data[column][i] = (*col_depth)[i];
      props.depth_c.data[column][i] = (*col_depth_c)[i];
      props.dz.data[column][i] = (*col_dz)[i];

      if (has_gas) {
        props.gas_saturation.data[column][i] = (*col_g_sat)[i];
        state.gas_density.data[column][i] = (*col_g_dens)[i];
      }

      if (has_ice) {
        state.ice_density.data[column][i] = (*col_i_dens)[i];
        props.ice_saturation.data[column][i] = (*col_i_sat)[i];
      }

      state.temperature.data[column][i] = (*col_temp)[i];
      props.thermal_conductivity.data[column][i] = (*col_cond)[i];
    }
    //fill surface variables

    state.surface_energy_source.data[column] = surface_energy_source[column];
    state.surface_water_source.data[column] = surface_water_source[column];
    state.snow_depth.data[column] = snow_depth[column];

    props.shortwave_radiation.data[column] = shortwave_radiation[column];
    props.longwave_radiation.data[column] = longwave_radiation[column];
    props.air_temperature.data[column] = air_temperature[column];
    props.vapor_pressure_air.data[column] = vapor_pressure_air[column];
    props.wind_speed.data[column] = wind_speed[column];
    props.precipitation.data[column] = precipitation[column];
    props.elevation.data[column] = elevation[column];
    props.aspect.data[column] = aspect[column];
    props.slope.data[column] = slope[column];

    for (int i = 0; i < state.total_component_concentration.columns; i++) {
      for (int j = 0; j < state.total_component_concentration.cells; j++) {
        for (int k = 0; k < state.total_component_concentration.components; k++) {
        }
      }
    }
  }

  //Fill the atmospheric abundances
  //NOTE: probably want to add an if statement here to only do this only once
  props.atm_n2 = atm_n2_;
  props.atm_o2 = atm_o2_;
  props.atm_co2 = atm_co2_;
  props.atm_ch4 = atm_ch4_;
  props.atm_n2o = atm_n2o_;
  props.atm_h2 = atm_h2_;
  props.atm_nh3 = atm_nh3_;
  props.heat_capacity = c_m_;
  props.field_capacity = pressure_at_field_capacity;
  props.wilting_point = pressure_at_wilting_point;

  Teuchos::OSTab tab = vo_->getOSTab();
  
  /*for (int column=0; column!=num_columns_local; ++column) {
  	*vo_->os() << "for column: " << column << std::endl;
	for (int i=0; i < ncells_per_col_; ++i) {
	   *vo_->os() << "T["<< i << "] = " << state.temperature.data[column][i] << std::endl; 
    }	
  }*/
}

void EcoSIM::CopyFromEcoSIM_process(const int column,
                                   const BGCProperties& props,
                                   const BGCState& state,
                                   const BGCAuxiliaryData& aux_data,
                                  const Tag& water_tag)
{

  Epetra_MultiVector& tcc= *(S_->GetPtrW<CompositeVector>(tcc_key_, Tags::DEFAULT, "state")->ViewComponent("cell",false));
  int tcc_num = tcc.NumVectors();

  auto& porosity = *(*S_->GetW<CompositeVector>(porosity_key_, Tags::DEFAULT, porosity_key_).ViewComponent("cell",false))(0);
  auto& liquid_saturation = *(*S_->GetW<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT, saturation_liquid_key_).ViewComponent("cell",false))(0);
  auto& water_content = *(*S_->GetW<CompositeVector>(water_content_key_, Tags::DEFAULT, water_content_key_).ViewComponent("cell",false))(0);
  //auto& suction_head = *(*S_->GetW<CompositeVector>(suc_key_, Tags::DEFAULT, suc_key_).ViewComponent("cell",false))(0);
  auto& relative_permeability = *(*S_->GetW<CompositeVector>(relative_permeability_key_, Tags::DEFAULT, relative_permeability_key_).ViewComponent("cell",false))(0);
  auto& liquid_density = *(*S_->GetW<CompositeVector>(liquid_density_key_, Tags::DEFAULT, liquid_density_key_).ViewComponent("cell",false))(0);
  auto& rock_density = *(*S_->GetW<CompositeVector>(rock_density_key_, Tags::DEFAULT, rock_density_key_).ViewComponent("cell",false))(0);
  auto& cell_volume = *(*S_->GetW<CompositeVector>(cell_volume_key_, Tags::DEFAULT, cell_volume_key_).ViewComponent("cell",false))(0);
  auto& hydraulic_conductivity = *(*S_->GetW<CompositeVector>(hydraulic_conductivity_key_, Tags::DEFAULT, hydraulic_conductivity_key_).ViewComponent("cell",false))(0);
  auto& bulk_density = *(*S_->GetW<CompositeVector>(bulk_density_key_, Tags::DEFAULT, bulk_density_key_).ViewComponent("cell",false))(0);

  //auto& surface_energy_source = *(*S_->GetW<CompositeVector>(surface_energy_source_key_, Tags::DEFAULT, surface_energy_source_key_).ViewComponent("cell", false))(0);
  auto& surface_energy_source = *(*S_->GetW<CompositeVector>(surface_energy_source_ecosim_key_, Tags::DEFAULT, surface_energy_source_ecosim_key_).ViewComponent("cell", false))(0);
  auto& subsurface_energy_source = *(*S_->GetW<CompositeVector>(subsurface_energy_source_key_, Tags::DEFAULT, subsurface_energy_source_key_).ViewComponent("cell", false))(0);

  //auto& surface_water_source = *(*S_->GetW<CompositeVector>(surface_water_source_key_, Tags::DEFAULT, surface_water_source_key_).ViewComponent("cell", false))(0);
  auto& surface_water_source = *(*S_->GetW<CompositeVector>(surface_water_source_ecosim_key_, Tags::DEFAULT, surface_water_source_ecosim_key_).ViewComponent("cell", false))(0);
  auto& subsurface_water_source = *(*S_->GetW<CompositeVector>(subsurface_water_source_key_, Tags::DEFAULT, subsurface_water_source_key_).ViewComponent("cell", false))(0);

  //auto& surface_test = *(*S_->GetW<CompositeVector>(surface_test_key_, Tags::DEFAULT, surface_test_key_).ViewComponent("cell", false))(0);
  auto& snow_depth = *(*S_->GetW<CompositeVector>(snow_depth_key_, Tags::DEFAULT, snow_depth_key_).ViewComponent("cell", false))(0);

  auto col_porosity = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_suc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_relative_permeability = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_b_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_energy_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_ss_water_source = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));

  //Gather columns on this process:
  num_columns_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  num_columns_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  //Trying to loop over processors now:
  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  num_columns_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  double energy_source_tot = state.surface_energy_source.data[0];
  double water_source_tot = state.surface_water_source.data[0];
  double snow_depth_cell = state.snow_depth.data[0];

  Teuchos::OSTab tab = vo_->getOSTab();	
  *vo_->os() << "snow depth (ATS): " << snow_depth_cell << " m" <<std::endl;

  //Loop over columns on this process
  for (int col=0; col!=num_columns_local; ++col) {

    if (has_gas) {
      auto& gas_saturation = *(*S_->GetW<CompositeVector>(saturation_gas_key_, Tags::DEFAULT, saturation_gas_key_).ViewComponent("cell",false))(0);
      auto& gas_density = *(*S_->GetW<CompositeVector>(gas_density_key_, Tags::DEFAULT, gas_density_key_).ViewComponent("cell",false))(0);
      for (int i=0; i < ncells_per_col_; ++i) {
        (*col_g_dens)[i] = state.gas_density.data[column][i];
        (*col_g_sat)[i] = props.gas_saturation.data[column][i];
      }

      ColumnToField_(column,gas_saturation,col_g_sat.ptr());
      ColumnToField_(column,gas_density,col_g_dens.ptr());
    }

    if (has_ice) {
      auto& ice_saturation = *(*S_->GetW<CompositeVector>(saturation_ice_key_, Tags::DEFAULT, saturation_ice_key_).ViewComponent("cell",false))(0);
      auto& ice_density = *(*S_->GetW<CompositeVector>(ice_density_key_, Tags::DEFAULT, ice_density_key_).ViewComponent("cell",false))(0);

      for (int i=0; i < ncells_per_col_; ++i) {
        (*col_i_dens)[i] = state.ice_density.data[column][i];
        (*col_i_sat)[i] = props.ice_saturation.data[column][i];
      }

      ColumnToField_(column,ice_saturation,col_i_sat.ptr());
      ColumnToField_(column,ice_density,col_i_dens.ptr());
    }

    auto& temp = *(*S_->GetW<CompositeVector>(T_key_, Tags::DEFAULT, "subsurface energy").ViewComponent("cell",false))(0);
    auto& thermal_conductivity = *(*S_->GetW<CompositeVector>(thermal_conductivity_key_, Tags::DEFAULT, thermal_conductivity_key_).ViewComponent("cell",false))(0);

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_temp)[i] = state.temperature.data[column][i];
      (*col_cond)[i] = props.thermal_conductivity.data[column][i];
    }

    ColumnToField_(column,temp, col_temp.ptr());
    ColumnToField_(column,thermal_conductivity,col_cond.ptr());

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_l_dens)[i] = state.liquid_density.data[column][i];
      (*col_porosity)[i] = state.porosity.data[column][i];
      (*col_wc)[i] = state.water_content.data[column][i];
      (*col_h_cond)[i] = state.hydraulic_conductivity.data[column][i];
      (*col_b_dens)[i] = state.bulk_density.data[column][i];

      (*col_ss_water_source)[i] = state.subsurface_water_source.data[column][i];
      (*col_ss_energy_source)[i] = state.subsurface_energy_source.data[column][i];
    }

    //double energy_source_tot = state.surface_energy_source.data[column];

    //As EcoSIM is hourly, but ATS is per second we need to divide the source by seconds per hour
    surface_energy_source[column] = state.surface_energy_source.data[column]/(3600.0);
    surface_water_source[column] = state.surface_water_source.data[column]/(3600.0);
    //surface_test[column] = state.surface_test.data[column];
    snow_depth[column] = state.snow_depth.data[column];

    //AG - 7/1/24 I'm not sure what this is doing I think it's from when I was testing comparing old to new values 
    //auto& new_e_source = *(*S_->GetW<CompositeVector>(surface_energy_source_key_, Tags::DEFAULT, surface_energy_source_key_).ViewComponent("cell", false))(0);
    //auto& new_w_source = *(*S_->GetW<CompositeVector>(surface_water_source_key_, Tags::DEFAULT, surface_water_source_key_).ViewComponent("cell", false))(0);

    ColumnToField_(column,liquid_saturation,col_l_sat.ptr());
    ColumnToField_(column,water_content,col_wc.ptr());
    ColumnToField_(column,relative_permeability,col_relative_permeability.ptr());
    ColumnToField_(column,hydraulic_conductivity,col_h_cond.ptr());
    ColumnToField_(column,bulk_density,col_b_dens.ptr());
    //ColumnToField_(column,plant_wilting_factor,col_wp.ptr());
    //ColumnToField_(column,rooting_depth_fraction,col_rf.ptr());
  }
}

/*
int EcoSIM::InitializeSingleColumn(int col)
{
  CopyToEcoSIM(column, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  //ecosim_datatest_wrapper(column, &bgc_props_, &bgc_sizes_);
  bgc_engine_->DataTest();

  int num_iterations = 1;

  bgc_engine_->Setup(bgc_props_, bgc_state_, bgc_sizes_, num_iterations, col);
  CopyEcoSIMStateToAmanzi(column, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

}

int EcoSIM::AdvanceSingleColumn(double dt, int col)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM(column, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 1;

 bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_sizes_, num_iterations, col);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyEcoSIMStateToAmanzi(column,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
 }
*/
int EcoSIM::InitializeSingleProcess(int proc)
{
  int num_iterations = 1;
  int num_columns = 1;
  
  num_columns = num_columns_local;

  Teuchos::OSTab tab = vo_->getOSTab();

  CopyToEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  //ecosim_datatest_wrapper(column, &bgc_props_, &bgc_sizes_);
  //bgc_engine_->DataTest();

  /*need some sort of assertions here to double check that the data is actually
  What I want it to be*/

  //Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "num_columns: " << num_columns << std::endl;
  *vo_->os() << "ncells_per_col_: " << ncells_per_col_ << std::endl;

  bgc_engine_->Setup(bgc_props_, bgc_state_, bgc_sizes_, num_iterations, num_columns,ncells_per_col_);
  CopyFromEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);
}

int EcoSIM::AdvanceSingleProcess(double dt, int proc)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 1;
  int num_columns = 1;

  num_columns = num_columns_local;

  bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_sizes_, num_iterations, num_columns);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyFromEcoSIM_process(proc,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
}

} // namespace EcoSIM
} // namespace Amanzi
