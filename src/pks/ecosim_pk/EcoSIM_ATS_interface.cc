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
    domain_surf_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

    // transport
    tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");
    //Remember tcc components are accessed by tcc[i][c] where i is the component and c is the cell

    //Flow
    poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_liquid_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    saturation_gas_key_ = Keys::readKey(*plist_,domain_,"saturation gas", "saturation_gas");
    saturation_ice_key_ = Keys::readKey(*plist_,domain_,"saturation ice", "saturation_ice");
    water_content_key_ = Keys::readKey(*plist_,domain_,"water content","water_content");
    rel_perm_key_ = Keys::readKey(*plist_,domain_,"relative permeability","relative_permeability");
    suc_key_ = Keys::readKey(*plist_,domain_,"suction","suction_head");
    liquid_den_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
    ice_den_key_ = Keys::readKey(*plist_, domain_, "mass density ice", "mass_density_ice");
    gas_den_key_ = Keys::readKey(*plist_, domain_,"mass density gas", "mass_density_gas");
    gas_den_key_test_ = Keys::readKey(*plist_, domain_, "mass density gas", "mass_density_gas");
    rock_den_key_ = Keys::readKey(*plist_, domain_, "density rock", "density_rock");

    //energy
    T_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
    therm_cond_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");

    //Other
    cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
    min_vol_frac_key_ = Keys::readKey(*plist_, domain_, "mineral volume fractions", "mineral_volume_fractions");
    ecosim_aux_data_key_ = Keys::readKey(*plist_, domain_, "ecosim aux data", "ecosim_aux_data");
    f_wp_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    f_root_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");

    //Custom Evaluator keys
    hydra_cond_key_ = Keys::readKey(*plist_, domain_, "hydraulic conductivity", "hydraulic_conductivity");
    bulk_dens_key_ = Keys::readKey(*plist_, domain_, "bulk density", "bulk_density");

    //Surface balance items
    sw_key_ =
      Keys::readKey(*plist_, domain_surf_, "incoming shortwave radiation", "incoming_shortwave_radiation");
    lw_key_ =
      Keys::readKey(*plist_,domain_surf_, "incoming longwave radiation", "incoming_longwave_radiation");
    air_temp_key_ = Keys::readKey(*plist_, domain_surf_, "air temperature", "air_temperature");
    vp_air_key_ = Keys::readKey(*plist_, domain_surf_, "vapor pressure air", "vapor_pressure_air");
    wind_speed_key_ = Keys::readKey(*plist_, domain_surf_, "wind speed", "wind_speed");
    prain_key_ = Keys::readKey(*plist_, domain_surf_, "precipitation rain", "precipitation_rain");
    elev_key_ = Keys::readKey(*plist_, domain_surf_, "elevation", "elevation");
    aspect_key_ = Keys::readKey(*plist_, domain_surf_, "aspect", "aspect");
    slope_key_ = Keys::readKey(*plist_, domain_surf_, "slope", "slope_magnitude");

    //Atmospheric abundance keys
    atm_n2_ = plist_->get<double>("atmospheric N2");
    atm_o2_ = plist_->get<double>("atmospheric O2");
    atm_co2_ = plist_->get<double>("atmospheric CO2");
    atm_ch4_ = plist_->get<double>("atmospheric CH4");
    atm_n2o_ = plist_->get<double>("atmospheric N2O");
    atm_h2_ = plist_->get<double>("atmospheric H2");
    atm_nh3_ = plist_->get<double>("atmospheric NH3");

    dt_ = plist_->get<double>("initial time step", 1.);
    c_m_ = plist_->get<double>("heat capacity [J mol^-1 K^-1]");

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
  mesh_surf_ = S_->GetMesh(domain_surf_);
  num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (unsigned int col = 0; col != num_cols_; ++col) {
    int f = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
    auto& col_iter = mesh_->cells_of_column(col);
    std::size_t ncol_cells = col_iter.size();

    // unclear which this should be:
    // -- col area is the true face area
    double col_area = mesh_->face_area(f);
    // -- col area is the projected face area
    // double col_area = mesh_surf_->cell_volume(col);

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
  requireAtNext(hydra_cond_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtCurrent(hydra_cond_key_, tag_current_, *S_, name_);

  requireAtNext(bulk_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtCurrent(bulk_dens_key_, tag_current_, *S_, name_);


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

  num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "columns on processor: " << num_cols_ << std::endl;

  //Now we call the engine's init state function which allocates the data
  bgc_engine_->InitState(bgc_props_, bgc_state_, bgc_aux_data_, ncells_per_col_, tcc_num, num_cols_);

  int ierr = 0;

  *vo_->os() << "printing bool:" << std::endl;
  *vo_->os() << (S_->HasRecord(suc_key_, Tags::DEFAULT)) << std::endl;
  if (S_->HasRecord(suc_key_, Tags::DEFAULT)) {
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
  }

  // Ensure dependencies are filled
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_wp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_root_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(suc_key_, Tags::DEFAULT).Update(*S_, name_);

  //Surface properties from met data
  S_->GetEvaluator(sw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(lw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(prain_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(air_temp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(vp_air_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(wind_speed_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(aspect_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(slope_key_, Tags::DEFAULT).Update(*S_, name_);

  if (S_->HasRecord(gas_den_key_test_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found mass density gas key" << std::endl;
    S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(gas_den_key_, Tags::DEFAULT).Update(*S_, name_);
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
    S_->GetEvaluator(therm_cond_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);
    has_energy = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Did not find temp key" << std::endl;
    has_energy = false;
  }

  if (S_->HasRecord(ice_den_key_, Tags::DEFAULT)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "found ice key" << std::endl;
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
    has_ice = true;
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Did not find ice key" << std::endl;
    has_ice = false;
  }

  //Initialize owned evaluators
  S_->GetW<CompositeVector>(hydra_cond_key_, Tags::DEFAULT, "hydraulic_conductivity").PutScalar(1.0);
  S_->GetRecordW(hydra_cond_key_, Tags::DEFAULT, "hydraulic_conductivity").set_initialized();

  S_->GetW<CompositeVector>(bulk_dens_key_, Tags::DEFAULT, "bulk_density").PutScalar(1.0);
  S_->GetRecordW(bulk_dens_key_, Tags::DEFAULT, "bulk_density").set_initialized();

  int num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //Looping over the columns and initializing
  /*for (int col=0; col!=num_cols_; ++col) {
    ierr = InitializeSingleColumn(col);
  }*/

  //loop over processes instead:
  ncols_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncols_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  *vo_->os() << "total columns cell_map: " << ncols_global << std::endl;
  *vo_->os() << "total columns from num_entities: " << ncols_global << std::endl;
  //Trying to loop over processors now:
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      std::cout << "on processor " << p_rank << std::endl;
      ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      std::cout << ncols_local << std::endl;

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
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_wp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(f_root_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(suc_key_, Tags::DEFAULT).Update(*S_, name_);

  //Surface data from met data
  S_->GetEvaluator(sw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(lw_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(air_temp_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(vp_air_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(wind_speed_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(prain_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(aspect_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(slope_key_, Tags::DEFAULT).Update(*S_, name_);

  if (has_gas) {
    S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(gas_den_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  if (has_ice) {
    S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  if (has_energy) {
    S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
    S_->GetEvaluator(therm_cond_key_, Tags::DEFAULT).Update(*S_, name_);
  }

  //Update owned evaluators
  Teuchos::RCP<const CompositeVector> hydra_cond = S_->GetPtr<CompositeVector>(hydra_cond_key_, Tags::DEFAULT);
  S_->GetEvaluator(hydra_cond_key_, Tags::DEFAULT).Update(*S_, name_);

  Teuchos::RCP<const CompositeVector> bulk_dens = S_->GetPtr<CompositeVector>(bulk_dens_key_, Tags::DEFAULT);
  S_->GetEvaluator(bulk_dens_key_, Tags::DEFAULT).Update(*S_, name_);

  AmanziMesh::Entity_ID num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

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

  /*S_->GetEvaluator("suction_head", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& suction_head = *(*S_->Get<CompositeVector>("suction_head", tag_next_)
          .ViewComponent("cell",false))(0);*/

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

  if (has_energy) {
    S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& temp = *(*S_->Get<CompositeVector>("temperature", tag_next_)
        .ViewComponent("cell",false))(0);

    S_->GetEvaluator("thermal_conductivity", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& thermal_conductivity = *(*S_->Get<CompositeVector>("thermal_conductivity", tag_next_)
            .ViewComponent("cell",false))(0);
  }

  //loop over processes instead:
  ncols_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncols_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  *vo_->os() << "total columns cell_map: " << ncols_global << std::endl;
  *vo_->os() << "total columns from num_entities: " << ncols_global << std::endl;
  //Trying to loop over processors now:
  int numProcesses, p_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  for (int k = 0; k < numProcesses; ++k) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (p_rank==k) {
      std::cout << "on processor " << p_rank << std::endl;
      ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      std::cout << ncols_local << std::endl;

      AdvanceSingleProcess(dt, p_rank);
    }
  }

}

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = vec[vec_index];
  }
}

void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Teuchos::Ptr<Epetra_SerialDenseVector> vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];

    (*col_vec)[i] = (*vec)[vec_index];
  }
}

void EcoSIM::MatrixFieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr)
  {
    int n_comp = m_arr.NumVectors();
    auto& col_iter = mesh_->cells_of_column(col);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        (*col_arr)(i,j) = m_arr[j][col_iter[i]];
      }
    }
  }

// helper function for pushing column back to field
void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID col, Epetra_Vector& vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    vec[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID col, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*vec)[col_iter[i]] = (*col_vec)[i];
  }
}

void EcoSIM::MatrixColumnToField_(AmanziMesh::Entity_ID col, Epetra_MultiVector& m_arr,
  Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr) {

    int n_comp = m_arr.NumVectors();
    auto& col_iter = mesh_->cells_of_column(col);

    for (int j=0; j!=n_comp; ++j){
      for (std::size_t i=0; i!=col_iter.size(); ++i) {
        m_arr[j][col_iter[i]] = (*col_arr)(i,j);
      }
    }

  }

// helper function for collecting column dz and depth
void EcoSIM::ColDepthDz_(AmanziMesh::Entity_ID col,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
  auto& col_iter = mesh_->cells_of_column(col);
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
  const Epetra_Vector& porosity = *(*S_->Get<CompositeVector>(poro_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_MultiVector& tcc= *(S_->GetPtr<CompositeVector>(tcc_key_, water_tag)->ViewComponent("cell"));
  int tcc_num = tcc.NumVectors();

  const Epetra_Vector& liquid_saturation = *(*S_->Get<CompositeVector>(saturation_liquid_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& water_content = *(*S_->Get<CompositeVector>(water_content_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& relative_permeability = *(*S_->Get<CompositeVector>(rel_perm_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_density = *(*S_->Get<CompositeVector>(liquid_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rock_density = *(*S_->Get<CompositeVector>(rock_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& cell_volume = *(*S_->Get<CompositeVector>(cv_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& hydraulic_conductivity = *(*S_->Get<CompositeVector>(hydra_cond_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& suction_head = *(*S_->Get<CompositeVector>(suc_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& bulk_density = *(*S_->Get<CompositeVector>(bulk_dens_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rooting_depth_fraction = *(*S_->Get<CompositeVector>(f_root_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& plant_wilting_factor = *(*S_->Get<CompositeVector>(f_wp_key_, water_tag).ViewComponent("cell", false))(0);

  const Epetra_Vector& shortwave_radiation = *(*S_->Get<CompositeVector>(sw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& longwave_radiation = *(*S_->Get<CompositeVector>(lw_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& air_temperature = *(*S_->Get<CompositeVector>(air_temp_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& vapor_pressure_air = *(*S_->Get<CompositeVector>(vp_air_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& wind_speed= *(*S_->Get<CompositeVector>(wind_speed_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& precipitation = *(*S_->Get<CompositeVector>(prain_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& elevation = *(*S_->Get<CompositeVector>(elev_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& aspect = *(*S_->Get<CompositeVector>(aspect_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& slope = *(*S_->Get<CompositeVector>(slope_key_, water_tag).ViewComponent("cell", false))(0);

  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_suc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
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

  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));

  //Gather columns on this process:
  ncols_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncols_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  *vo_->os() << "total columns cell_map: " << ncols_global << std::endl;
  *vo_->os() << "total columns from num_entities: " << ncols_global_ptype << std::endl;
  //Trying to loop over processors now:
  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "on processor " << p_rank << std::endl;
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  std::cout << ncols_local << std::endl;

  //Loop over columns on this process
  for (int col=0; col!=ncols_local; ++col) {
    FieldToColumn_(col,porosity,col_poro.ptr());
    FieldToColumn_(col,liquid_saturation,col_l_sat.ptr());
    FieldToColumn_(col,water_content,col_wc.ptr());
    FieldToColumn_(col,relative_permeability,col_rel_perm.ptr());
    FieldToColumn_(col,liquid_density,col_l_dens.ptr());
    FieldToColumn_(col,rock_density,col_r_dens.ptr());
    FieldToColumn_(col,cell_volume,col_vol.ptr());
    FieldToColumn_(col,hydraulic_conductivity,col_h_cond.ptr());
    FieldToColumn_(col,bulk_density,col_b_dens.ptr());
    FieldToColumn_(col,plant_wilting_factor,col_wp.ptr());
    FieldToColumn_(col,rooting_depth_fraction,col_rf.ptr());

    MatrixFieldToColumn_(col, tcc, col_tcc.ptr());

    if (has_gas) {
      const Epetra_Vector& gas_saturation = *(*S_->Get<CompositeVector>(saturation_gas_key_, water_tag).ViewComponent("cell", false))(0);
      const Epetra_Vector& gas_density = *(*S_->Get<CompositeVector>(gas_den_key_, water_tag).ViewComponent("cell", false))(0);

      FieldToColumn_(col,gas_saturation,col_g_sat.ptr());
      FieldToColumn_(col,gas_density,col_g_dens.ptr());
    }

    if (has_ice) {
      const Epetra_Vector& ice_saturation = *(*S_->Get<CompositeVector>(saturation_ice_key_, water_tag).ViewComponent("cell", false))(0);
      const Epetra_Vector& ice_density = *(*S_->Get<CompositeVector>(ice_den_key_, water_tag).ViewComponent("cell", false))(0);

      FieldToColumn_(col,ice_saturation,col_i_sat.ptr());
      FieldToColumn_(col,ice_density,col_i_dens.ptr());
    }

    if (has_energy) {
      const Epetra_Vector& temp = *(*S_->Get<CompositeVector>(T_key_, water_tag).ViewComponent("cell", false))(0);
      const Epetra_Vector& thermal_conductivity = *(*S_->Get<CompositeVector>(therm_cond_key_, water_tag).ViewComponent("cell", false))(0);

      FieldToColumn_(col,temp, col_temp.ptr());
      FieldToColumn_(col,thermal_conductivity,col_cond.ptr());
    }

    // This is for computing depth
    ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

    for (int i=0; i < ncells_per_col_; ++i) {
      state.liquid_density.data[col][i] = (*col_l_dens)[i];
      state.porosity.data[col][i] = (*col_poro)[i];
      state.water_content.data[col][i] = (*col_wc)[i];
      state.hydraulic_conductivity.data[col][i] = (*col_h_cond)[i];
      state.bulk_density.data[col][i] = (*col_b_dens)[i];
      //state.suction_head.data[i] = (*col_suc)[i];
      props.plant_wilting_factor.data[col][i] = (*col_wp)[i];
      props.rooting_depth_fraction.data[col][i] = (*col_rf)[i];
      props.liquid_saturation.data[col][i] = (*col_l_sat)[i];
      props.relative_permeability.data[col][i] = (*col_rel_perm)[i];
      props.volume.data[col][i] = (*col_vol)[i];
      props.depth.data[col][i] = (*col_depth)[i];
      props.dz.data[col][i] = (*col_dz)[i];

      if (has_gas) {
        props.gas_saturation.data[col][i] = (*col_g_sat)[i];
        state.gas_density.data[col][i] = (*col_g_dens)[i];
      }

      if (has_ice) {
        state.ice_density.data[col][i] = (*col_i_dens)[i];
        props.ice_saturation.data[col][i] = (*col_i_sat)[i];
      }

      if (has_energy) {
        state.temperature.data[col][i] = (*col_temp)[i];
        props.thermal_conductivity.data[col][i] = (*col_cond)[i];
      }
    }

    *vo_->os() << "filling surface props" << std::endl;
    *vo_->os() << "size of shortwave in struct is: " << props.shortwave_radiation.size << std::endl; 
    *vo_->os() << "size of shortwave in state is: " << shortwave_radiation.MyLength() << std::endl;

    //fill surface variables
    props.shortwave_radiation.data[col] = shortwave_radiation[col];
    props.longwave_radiation.data[col] = longwave_radiation[col];
    props.air_temperature.data[col] = air_temperature[col];
    props.vapor_pressure_air.data[col] = vapor_pressure_air[col];
    props.wind_speed.data[col] = wind_speed[col];
    props.precipitation.data[col] = precipitation[col];
    props.elevation.data[col] = elevation[col];
    props.aspect.data[col] = aspect[col];
    props.slope.data[col] = slope[col];

    for (int component=0; component < tcc_num; ++component) {
      for (int i=0; i < ncells_per_col_; ++i) {
        state.total_component_concentration.data[col][component][i] = (*col_tcc)(i,component);
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

}

void EcoSIM::CopyFromEcoSIM_process(const int col,
                                   const BGCProperties& props,
                                   const BGCState& state,
                                   const BGCAuxiliaryData& aux_data,
                                  const Tag& water_tag)
{

  Epetra_MultiVector& tcc= *(S_->GetPtrW<CompositeVector>(tcc_key_, Amanzi::Tags::NEXT, "state")->ViewComponent("cell",false));
  int tcc_num = tcc.NumVectors();

  auto& porosity = *(*S_->GetW<CompositeVector>(poro_key_, Amanzi::Tags::NEXT, poro_key_).ViewComponent("cell",false))(0);
  auto& liquid_saturation = *(*S_->GetW<CompositeVector>(saturation_liquid_key_, Amanzi::Tags::NEXT, saturation_liquid_key_).ViewComponent("cell",false))(0);
  auto& water_content = *(*S_->GetW<CompositeVector>(water_content_key_, Amanzi::Tags::NEXT, water_content_key_).ViewComponent("cell",false))(0);
  //auto& suction_head = *(*S_->GetW<CompositeVector>(suc_key_, Amanzi::Tags::NEXT, suc_key_).ViewComponent("cell",false))(0);
  auto& relative_permeability = *(*S_->GetW<CompositeVector>(rel_perm_key_, Amanzi::Tags::NEXT, rel_perm_key_).ViewComponent("cell",false))(0);
  auto& liquid_density = *(*S_->GetW<CompositeVector>(liquid_den_key_, Amanzi::Tags::NEXT, liquid_den_key_).ViewComponent("cell",false))(0);
  auto& rock_density = *(*S_->GetW<CompositeVector>(rock_den_key_, Amanzi::Tags::NEXT, rock_den_key_).ViewComponent("cell",false))(0);
  auto& cell_volume = *(*S_->GetW<CompositeVector>(cv_key_, Amanzi::Tags::NEXT, cv_key_).ViewComponent("cell",false))(0);
  auto& hydraulic_conductivity = *(*S_->GetW<CompositeVector>(hydra_cond_key_, Amanzi::Tags::NEXT, hydra_cond_key_).ViewComponent("cell",false))(0);
  auto& bulk_density = *(*S_->GetW<CompositeVector>(bulk_dens_key_, Amanzi::Tags::NEXT, bulk_dens_key_).ViewComponent("cell",false))(0);

  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_suc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
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

  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));

  //Gather columns on this process:
  ncols_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncols_global_ptype = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  *vo_->os() << "total columns cell_map: " << ncols_global << std::endl;
  *vo_->os() << "total columns from num_entities: " << ncols_global_ptype << std::endl;
  //Trying to loop over processors now:
  int p_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "on processor " << p_rank << std::endl;
  ncols_local = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  std::cout << ncols_local << std::endl;

  //Loop over columns on this process
  for (int col=0; col!=ncols_local; ++col) {
    if (has_gas) {
      auto& gas_saturation = *(*S_->GetW<CompositeVector>(saturation_gas_key_, Amanzi::Tags::NEXT, saturation_gas_key_).ViewComponent("cell",false))(0);
      auto& gas_density = *(*S_->GetW<CompositeVector>(gas_den_key_, Amanzi::Tags::NEXT, gas_den_key_).ViewComponent("cell",false))(0);
      for (int i=0; i < ncells_per_col_; ++i) {
        (*col_g_dens)[i] = state.gas_density.data[col][i];
        (*col_g_sat)[i] = props.gas_saturation.data[col][i];
      }

      FieldToColumn_(col,gas_saturation,col_g_sat.ptr());
      FieldToColumn_(col,gas_density,col_g_dens.ptr());
    }

    if (has_ice) {
      auto& ice_saturation = *(*S_->GetW<CompositeVector>(saturation_ice_key_, Amanzi::Tags::NEXT, saturation_ice_key_).ViewComponent("cell",false))(0);
      auto& ice_density = *(*S_->GetW<CompositeVector>(ice_den_key_, Amanzi::Tags::NEXT, ice_den_key_).ViewComponent("cell",false))(0);

      for (int i=0; i < ncells_per_col_; ++i) {
        (*col_i_dens)[i] = state.ice_density.data[col][i];
        (*col_i_sat)[i] = props.ice_saturation.data[col][i];
      }

      ColumnToField_(col,ice_saturation,col_i_sat.ptr());
      ColumnToField_(col,ice_density,col_i_dens.ptr());
    }

    if (has_energy) {
      auto& temp = *(*S_->GetW<CompositeVector>(T_key_, Amanzi::Tags::NEXT, "subsurface energy").ViewComponent("cell",false))(0);
      auto& thermal_conductivity = *(*S_->GetW<CompositeVector>(therm_cond_key_, Amanzi::Tags::NEXT, therm_cond_key_).ViewComponent("cell",false))(0);

      for (int i=0; i < ncells_per_col_; ++i) {
        (*col_temp)[i] = state.temperature.data[col][i];
        (*col_cond)[i] = props.thermal_conductivity.data[col][i];
      }

      ColumnToField_(col,temp, col_temp.ptr());
      ColumnToField_(col,thermal_conductivity,col_cond.ptr());
    }

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_l_dens)[i] = state.liquid_density.data[col][i];
      (*col_poro)[i] = state.porosity.data[col][i];
      (*col_wc)[i] = state.water_content.data[col][i];
      (*col_h_cond)[i] = state.hydraulic_conductivity.data[col][i];
      (*col_b_dens)[i] = state.bulk_density.data[col][i];

      if (has_gas) {
        (*col_g_dens)[i] = state.gas_density.data[col][i];
      }

      if (has_ice) {
        (*col_i_dens)[i] = state.ice_density.data[col][i];
      }

      if (has_energy) {
        (*col_temp)[i] = state.temperature.data[col][i];
      }
    }

    ColumnToField_(col,porosity,col_poro.ptr());
    ColumnToField_(col,liquid_saturation,col_l_sat.ptr());
    ColumnToField_(col,water_content,col_wc.ptr());
    ColumnToField_(col,relative_permeability,col_rel_perm.ptr());
    ColumnToField_(col,liquid_density,col_l_dens.ptr());
    ColumnToField_(col,rock_density,col_r_dens.ptr());
    ColumnToField_(col,cell_volume,col_vol.ptr());
    ColumnToField_(col,hydraulic_conductivity,col_h_cond.ptr());
    ColumnToField_(col,bulk_density,col_b_dens.ptr());
    //ColumnToField_(col,plant_wilting_factor,col_wp.ptr());
    //ColumnToField_(col,rooting_depth_fraction,col_rf.ptr());
  }
}
/*
int EcoSIM::InitializeSingleColumn(int col)
{
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  //ecosim_datatest_wrapper(col, &bgc_props_, &bgc_sizes_);
  bgc_engine_->DataTest();

  int num_iterations = 1;

  bgc_engine_->Setup(bgc_props_, bgc_state_, bgc_sizes_, num_iterations, col);
  CopyEcoSIMStateToAmanzi(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

}

int EcoSIM::AdvanceSingleColumn(double dt, int col)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 1;

 bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_sizes_, num_iterations, col);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyEcoSIMStateToAmanzi(col,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
 }
*/
int EcoSIM::InitializeSingleProcess(int proc)
{
  CopyToEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  //ecosim_datatest_wrapper(col, &bgc_props_, &bgc_sizes_);
  bgc_engine_->DataTest();

  int num_iterations = 1;
  int ncols = 1;

  bgc_engine_->Setup(bgc_props_, bgc_state_, bgc_sizes_, num_iterations, ncols);
  CopyFromEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

}

int EcoSIM::AdvanceSingleProcess(double dt, int proc)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM_process(proc, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 1;
  int ncols = 1;

  bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_sizes_, num_iterations, ncols);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyFromEcoSIM_process(proc,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
}

} // namespace EcoSIM
} // namespace Amanzi
