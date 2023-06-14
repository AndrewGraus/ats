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
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "Mesh.hh"

// include evaluators here
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
    //grab the relevant domains, surface is needed to find the columns later
    domain_ = plist_->get<std::string>("domain name", "domain");
    domain_surf_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

    // obtain key of fields
    // grid position (X,Y,Z) - can we just pass Z and assume X and Y are 0?
    // Aspect in geometric format
    // water table depth - There's a water table evaluator in
    // /src/constitutive_relations/column_integrators/ but I don't see it used
    // anywhere can we just use it here?


    //For now a few of the evaluators don't work because they are no initialized
    // I've commented out elevation and replaced others with something that does
    //exist:
    //
    // mass_density_ice
    // mass_gas_density
    //
    // additionally temperature doesn't work because it is owned by energy
    //
    //
    // Simple tests with the keys

    // transport
    tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");
    //Remember tcc components are accessed by tcc[i][c] where i is the component and c is the cell

    //Flow
    poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_liquid_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    saturation_gas_key_ = Keys::readKey(*plist_,domain_,"saturation gas", "saturation_gas");
    saturation_ice_key_ = Keys::readKey(*plist_,domain_,"saturation ice", "saturation_ice");
    //elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");
    water_content_key_ = Keys::readKey(*plist_,domain_,"water content","water_content");
    rel_perm_key_ = Keys::readKey(*plist_,domain_,"relative permeability","relative_permeability");

    //densities
    //If we need bulk density do we need volume fractions of each quantity?
    //This can be computed from the saturations and porosity (I think) via:
    // f_rock = (1 - porosity)
    // f_liq = S_liq * porosity
    // f_gas = S_gas * porosity
    // f_ice = S_ice * porosity

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

    //Evaluator keys
    hydra_cond_key_ = Keys::readKey(*plist_, domain_, "hydraulic conductivity", "hydraulic_conductivity");
    //bulk_dens_key_ = Keys::readKey(*plist_, domain_, "bulk density", "bulk_density");

    // parameters
    // initial timestep
    dt_ = plist_->get<double>("initial time step", 1.);
    //Heat capacity looks like the default units are molar heat capacity
    c_m_ = plist_->get<double>("heat capacity [J mol^-1 K^-1]");

    //They also sometimes use a version of heat capacity that is just this
    //quantity times 1e-6:
    //ka_ = 1.e-6 * plist_.get<double>("heat capacity [J kg^-1 K^-1]");
    //Unclear what we need

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

  /*requireAtNext(bulk_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtCurrent(bulk_dens_key_, tag_current_, *S_, name_);*/


  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Setup of PK was successful"
    << vo_->reset() << std::endl << std::endl;
  }
}

// -- Initialize owned (dependent) variables.
void EcoSIM::Initialize() {
  //Now we have to initalize the variables (i.e. give them initial values)
  //In our PK it will only be done for variables owned by the PK
  //Keeping an example of how it's done generically here:
  //S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_).PutScalar(0.);
  //S_->GetRecordW("co2_decomposition", tag_next_, name_).set_initialized();

  //In alquimia they initialize the axuiliary data via a function called InitializeCVField
  //which can be found in Amanzi/src/PKs/PK_Physical.cc:

  // initialize fields as soon as possible
  /*for (size_t i = 0; i < aux_names_.size(); ++i) {
    InitializeCVField(S_, *vo_, aux_names_[i], tag_next_, passwd_, 0.0);
  }*/

  //Need to know the number of components to initialize data structures
  const Epetra_MultiVector& tcc= *(S_->GetPtr<CompositeVector>(tcc_key_, Tags::DEFAULT)->ViewComponent("cell"));
  int tcc_num = tcc.NumVectors();

  //Now we call the engine's init state function which allocates the data
  bgc_engine_->InitState(bgc_props_, bgc_state_, bgc_aux_data_, ncells_per_col_, tcc_num);

  int ierr = 0;

  // Ensure dependencies are filled
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "testing keys" << std::endl;

  //Here we put the checks for the optional keys
  //Temperature, ice and gas
  //plist_->print(std::cout);

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
  //*vo_->os() << "Getting hydraulic conductivity" << std::endl;
  S_->GetW<CompositeVector>(hydra_cond_key_, Tags::DEFAULT, "hydraulic_conductivity").PutScalar(1.0);
  //*vo_->os() << "recording to hydraulic" << std::endl;
  S_->GetRecordW(hydra_cond_key_, Tags::DEFAULT, "hydraulic_conductivity").set_initialized();

  //S_->GetW<CompositeVector>(bulk_dens_key_, tag_next_, name_).PutScalar(1.0);
  //S_->GetRecordW(bulk_dens_key_, tag_next_, name_).set_initialized();

  int num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //Looping over the columns and initializing
  for (int col=0; col!=num_cols_; ++col) {
    ierr = InitializeSingleColumn(col);
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

  std::cout << "\nBegin Advance\n";
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
  //S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);

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

  //Teuchos::RCP<const CompositeVector> bulk_dens = S_->GetPtr<CompositeVector>(bulk_dens_key_, Tags::DEFAULT);
  //S_->GetEvaluator(bulk_dens_key_, Tags::DEFAULT).Update(*S_, name_);

  AmanziMesh::Entity_ID num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // grab the required fields
  S_->GetEvaluator("porosity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& porosity = *(*S_->Get<CompositeVector>("porosity", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_saturation = *(*S_->Get<CompositeVector>("saturation_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  //S_->GetEvaluator("elevation", tag_next_).Update(*S_, name_);
  //const Epetra_MultiVector& elevation = *S_->Get<CompositeVector>("elevation", tag_next_)
  //        .ViewComponent("cell",false);

  S_->GetEvaluator("water_content", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& water_content = *(*S_->Get<CompositeVector>("water_content", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("relative_permeability", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& relative_permeability = *(*S_->Get<CompositeVector>("relative_permeability", tag_next_)
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

  if (has_gas) {
    S_->GetEvaluator("mass_density_gas", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& gas_density = *(*S_->Get<CompositeVector>("mass_density_gas", tag_next_)
            .ViewComponent("cell",false))(0);

    S_->GetEvaluator("saturation_gas", tag_next_).Update(*S_, name_);
    const Epetra_MultiVector& gas_saturation = *(*S_->Get<CompositeVector>("saturation_gas", tag_next_)
            .ViewComponent("cell",false))(0);
  }

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

  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col=0; col!=num_cols_; ++col) {

    auto& col_iter = mesh_->cells_of_column(col);
    ncells_per_col_ = col_iter.size();

    //Copy to EcoSIM structures

    std::cout << "\nAdvancing col "<< col <<"\n";
    AdvanceSingleColumn(dt, col);
    std::cout << "\nfinished advancing column\n";

  } // end loop over columns

  // mark primaries as changed

  // Compute the next time step.
  // * will we need to do this? *
  //ComputeNextTimeStep();

  //return failed;

  std::cout << "\nEnd Advance\n";
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
void EcoSIM::CopyToEcoSIM(int col,
                                 BGCProperties& props,
                                 BGCState& state,
                                 BGCAuxiliaryData& aux_data,
                               const Tag& water_tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "Starting copy step" << std::endl;
  //Fill state with ATS variables that are going to be changed by EcoSIM
  const Epetra_Vector& porosity = *(*S_->Get<CompositeVector>(poro_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_MultiVector& tcc= *(S_->GetPtr<CompositeVector>(tcc_key_, water_tag)->ViewComponent("cell"));
  int tcc_num = tcc.NumVectors();

  const Epetra_Vector& liquid_saturation = *(*S_->Get<CompositeVector>(saturation_liquid_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& elevation = *S_->Get<CompositeVector>(elev_key_, water_tag).ViewComponent("cell", false);
  const Epetra_Vector& water_content = *(*S_->Get<CompositeVector>(water_content_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& relative_permeability = *(*S_->Get<CompositeVector>(rel_perm_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_density = *(*S_->Get<CompositeVector>(liquid_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rock_density = *(*S_->Get<CompositeVector>(rock_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& cell_volume = *(*S_->Get<CompositeVector>(cv_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& hydraulic_conductivity = *(*S_->Get<CompositeVector>(cv_key_, water_tag).ViewComponent("cell", false))(0);

  //Define the column vectors to hold the data
  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_elev = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));\
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));


  //For the concentration I do not want a vector but a matrix
  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));

  //Here is where we should do the various field-to-column calls to then pass along
  //to the data structures that will pass the data to EcoSIM
  //Format is:
  //FieldToColumn_(column index, dataset to copy from, vector to put the data in)
  FieldToColumn_(col,porosity,col_poro.ptr());
  FieldToColumn_(col,liquid_saturation,col_l_sat.ptr());
  //FieldToColumn_(col,elevation,col_elev.ptr());
  FieldToColumn_(col,water_content,col_wc.ptr());
  FieldToColumn_(col,relative_permeability,col_rel_perm.ptr());
  FieldToColumn_(col,liquid_density,col_l_dens.ptr());
  FieldToColumn_(col,rock_density,col_r_dens.ptr());
  FieldToColumn_(col,cell_volume,col_vol.ptr());
  FieldToColumn_(col,hydraulic_conductivity,col_h_cond.ptr());

  MatrixFieldToColumn_(col, tcc, col_tcc.ptr());

  /**vo_->os() << "Total Comp: " << tcc_num << std::endl;
  *vo_->os() << "Total cells: " << ncells_per_col_ << std::endl;
  for (int i=0; i < tcc_num; ++i) {
    Epetra_SerialDenseVector col_comp(ncells_per_col_);
    Epetra_SerialDenseVector tcc_comp(ncells_per_col_);
    *vo_->os() << "Component: " << i << std::endl;
    for (int j=0; j<ncells_per_col_; ++j){
      *vo_->os() << "Cell: " << j << std::endl;
      col_comp(j) = (*col_tcc)(i,j);
      tcc_comp[j] = tcc[i][j];
    }
    FieldToColumn_(col,Teuchos::ptr(&tcc_comp),Teuchos::ptr(&col_comp));
  }*/



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

  // I think I need to loop over the column data and save it to the data
  // structures. Eventually I could probably rewrite FieldToColumn_ to do this
  // have to fill tcc separately (I think)

  for (int j=0; j < tcc_num; ++j) {
    *vo_->os() << "component: "<< j << std::endl;
    for (int i=0; i < ncells_per_col_; ++i) {
      *vo_->os() << "cell: "<< i << std::endl;
      *vo_->os() << "col arr: "<< (*col_tcc)(i,j) << std::endl;
      *vo_->os() << "m_arr: "<< state.total_component_concentration.data[j][i] << std::endl;
      state.total_component_concentration.data[j][i] = (*col_tcc)(i,j);
    }
  }

  for (int i=0; i < ncells_per_col_; ++i) {
    state.liquid_density.data[i] = (*col_l_dens)[i];
    state.porosity.data[i] = (*col_poro)[i];
    state.water_content.data[i] = (*col_wc)[i];
    state.hydraulic_conductivity.data[i] = (*col_h_cond)[i];
    props.liquid_saturation.data[i] = (*col_l_sat)[i];
    //props.elevation.data[i] = (*col_elev)[i];
    props.relative_permeability.data[i] = (*col_rel_perm)[i];
    props.volume.data[i] = (*col_vol)[i];

    if (has_gas) {
      props.gas_saturation.data[i] = (*col_g_sat)[i];
      state.gas_density.data[i] = (*col_g_dens)[i];
    }

    if (has_ice) {
      state.ice_density.data[i] = (*col_i_dens)[i];
      props.ice_saturation.data[i] = (*col_i_sat)[i];
    }

    if (has_energy) {
      state.temperature.data[i] = (*col_temp)[i];
      props.thermal_conductivity.data[i] = (*col_cond)[i];
    }

  }
  //mat_props.volume = mesh_->cell_volume(cell;z
  //mat_props.saturation = water_saturation[0][cell];

  // Auxiliary data -- block copy.
  /*if (S_->HasRecord(bgc_aux_data_key_, tag_next_)) {
    aux_data_ = S_->GetW<CompositeVector>(bgc_aux_data_key_, tag_next_, passwd_).ViewComponent("cell");
    int num_aux_ints = bgc_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = bgc_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      aux_data.aux_ints.data[i] = (int)cell_aux_ints[cell];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      aux_data.aux_doubles.data[i] = cell_aux_doubles[cell];
    }
  }*/
}

void EcoSIM::CopyEcoSIMStateToAmanzi(
    const int col,
    const BGCProperties& props,
    const BGCState& state,
    const BGCAuxiliaryData& aux_data,
    const Tag& water_tag)
{
  CopyFromEcoSIM(col, props, state, aux_data, water_tag);
}

void EcoSIM::CopyFromEcoSIM(const int col,
                                   const BGCProperties& props,
                                   const BGCState& state,
                                   const BGCAuxiliaryData& aux_data,
                                  const Tag& water_tag)
{
  // If the chemistry has modified the porosity and/or density, it needs to
  // be updated here.
  // (this->water_density())[cell] = state.water_density;
  // (this->porosity())[cell] = state.porosity;

  Epetra_MultiVector& tcc= *(S_->GetPtrW<CompositeVector>(tcc_key_, Amanzi::Tags::NEXT, "state")->ViewComponent("cell",false));
  int tcc_num = tcc.NumVectors();

  auto& porosity = *(*S_->GetW<CompositeVector>(poro_key_, Amanzi::Tags::NEXT, poro_key_).ViewComponent("cell",false))(0);
  auto& liquid_saturation = *(*S_->GetW<CompositeVector>(saturation_liquid_key_, Amanzi::Tags::NEXT, saturation_liquid_key_).ViewComponent("cell",false))(0);
  //auto& elevation = S_->GetPtrW<CompositeVector>(elev_key_, Amanzi::Tags::NEXT, passwd_).ViewComponent("cell");
  auto& water_content = *(*S_->GetW<CompositeVector>(water_content_key_, Amanzi::Tags::NEXT, water_content_key_).ViewComponent("cell",false))(0);
  auto& relative_permeability = *(*S_->GetW<CompositeVector>(rel_perm_key_, Amanzi::Tags::NEXT, rel_perm_key_).ViewComponent("cell",false))(0);
  auto& liquid_density = *(*S_->GetW<CompositeVector>(liquid_den_key_, Amanzi::Tags::NEXT, liquid_den_key_).ViewComponent("cell",false))(0);
  auto& rock_density = *(*S_->GetW<CompositeVector>(rock_den_key_, Amanzi::Tags::NEXT, rock_den_key_).ViewComponent("cell",false))(0);
  auto& cell_volume = *(*S_->GetW<CompositeVector>(cv_key_, Amanzi::Tags::NEXT, cv_key_).ViewComponent("cell",false))(0);
  auto& hydraulic_conductivity = *(*S_->GetW<CompositeVector>(hydra_cond_key_, Amanzi::Tags::NEXT, hydra_cond_key_).ViewComponent("cell",false))(0);


  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_elev = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //For the concentration I do not want a vector but a matrix
  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseMatrix(tcc_num,ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_h_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  if (has_gas) {
    auto& gas_saturation = *(*S_->GetW<CompositeVector>(saturation_gas_key_, Amanzi::Tags::NEXT, saturation_gas_key_).ViewComponent("cell", false))(0);
    auto& gas_density = *(*S_->GetW<CompositeVector>(gas_den_key_, Amanzi::Tags::NEXT, gas_den_key_).ViewComponent("cell", false))(0);

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_g_dens)[i] = state.gas_density.data[i];
      (*col_g_sat)[i] = props.gas_saturation.data[i];
    }

    ColumnToField_(col,gas_saturation,col_g_sat.ptr());
    ColumnToField_(col,gas_density,col_g_dens.ptr());

  }

  if (has_ice) {
    auto& ice_saturation = *(*S_->GetW<CompositeVector>(saturation_ice_key_, Amanzi::Tags::NEXT, saturation_ice_key_).ViewComponent("cell",false))(0);
    auto& ice_density = *(*S_->GetW<CompositeVector>(ice_den_key_, Amanzi::Tags::NEXT, ice_den_key_).ViewComponent("cell",false))(0);

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_i_dens)[i] = state.ice_density.data[i];
      (*col_i_sat)[i] = props.ice_saturation.data[i];
    }

    ColumnToField_(col,ice_saturation,col_i_sat.ptr());
    ColumnToField_(col,ice_density,col_i_dens.ptr());
  }

  if (has_energy) {
    auto& temp = *(*S_->GetW<CompositeVector>(T_key_, Amanzi::Tags::NEXT, "energy").ViewComponent("cell",false))(0);
    auto& thermal_conductivity = *(*S_->GetW<CompositeVector>(therm_cond_key_, Amanzi::Tags::NEXT, therm_cond_key_).ViewComponent("cell",false))(0);

    for (int i=0; i < ncells_per_col_; ++i) {
      (*col_temp)[i] = state.temperature.data[i];
      (*col_cond)[i] = props.thermal_conductivity.data[i];
    }

    ColumnToField_(col,temp, col_temp.ptr());
    ColumnToField_(col,thermal_conductivity,col_cond.ptr());
  }

  for (int i=0; i < ncells_per_col_; ++i) {
    (*col_l_dens)[i] = state.liquid_density.data[i];
    (*col_poro)[i] = state.porosity.data[i];
    (*col_wc)[i] = state.water_content.data[i];
    (*col_l_sat)[i] = props.liquid_saturation.data[i];
    //(*col_elev)[i] = props.elevation.data[i];
    (*col_rel_perm)[i] = props.relative_permeability.data[i];
    (*col_cond)[i] = props.thermal_conductivity.data[i];
    (*col_h_cond)[i] = state.hydraulic_conductivity.data[i];
    (*col_vol)[i] = props.volume.data[i];
  }

  //Take new values from Ecosim state and put them into the secondary data structure
  //for backing back into amanzi state
  for (int j=0; j < tcc_num; ++j) {
    *vo_->os() << "component: "<< j << std::endl;
    for (int i=0; i < ncells_per_col_; ++i) {
      *vo_->os() << "cell: "<< i << std::endl;
      *vo_->os() << "col arr: "<< (*col_tcc)(i,j) << std::endl;
      *vo_->os() << "m_arr: "<< state.total_component_concentration.data[j][i] << std::endl;
      (*col_tcc)(i,j) = state.total_component_concentration.data[j][i];
    }
  }

  //Here is where the auxiliary data is filled need to try to change this to columns
  //This may not be trivial
  /*if (S_->HasRecord(bgc_aux_data_key_, tag_next_)) {
    aux_data_ = S_->GetW<CompositeVector>(bgc_aux_data_key_, tag_next_, passwd_).ViewComponent("cell");

    int num_aux_ints = bgc_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = bgc_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      cell_aux_ints[cell] = (double)aux_data.aux_ints.data[i];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      cell_aux_doubles[cell] = aux_data.aux_doubles.data[i];
    }
  }*/

  //pack this data back into the num_columns
  //ColumnToField_(col,tcc,col_tcc.ptr());
  ColumnToField_(col,porosity,col_poro.ptr());
  ColumnToField_(col,liquid_saturation,col_l_sat.ptr());
  //ColumnToField_(col,elevation,col_elev.ptr());
  ColumnToField_(col,water_content,col_wc.ptr());
  ColumnToField_(col,relative_permeability,col_rel_perm.ptr());
  ColumnToField_(col,liquid_density,col_l_dens.ptr());
  ColumnToField_(col,rock_density,col_r_dens.ptr());
  ColumnToField_(col,cell_volume,col_vol.ptr());
  ColumnToField_(col,hydraulic_conductivity,col_h_cond.ptr());

  MatrixColumnToField_(col, tcc, col_tcc.ptr());
}

/* *******************************************************************
* This helper performs initialization on a single column within Amanzi's state.
******************************************************************* */
int EcoSIM::InitializeSingleColumn(int col)
{
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  std::cout << "running data test" << std::endl;

  //ecosim_datatest_wrapper(col, &bgc_props_, &bgc_sizes_);
  bgc_engine_->DataTest();

  int num_iterations = 1;

  bgc_engine_->Setup(bgc_props_, bgc_state_, bgc_sizes_, num_iterations, col);
  CopyEcoSIMStateToAmanzi(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  // ETC: hacking to get consistent solution -- if there is no water
  // (e.g. surface system, we still need to call EnforceCondition() as it also
  // gets aux data set up correctly.  But the concentrations need to be
  // overwritten as 0 to get expected output.  Therefore we manually overwrite
  // this now.  Previously this happened due to a bug in ATS's reactive
  // transport coupler -- happy accidents.
  //if (alq_mat_props_.saturation <= saturation_tolerance_)
  //  for (int i=0; i!=aqueous_components_->NumVectors(); ++i) (*aqueous_components_)[i][cell] = 0.;
  //return 0;
}

/* *******************************************************************
* This helper advances the solution on a single cell within Amanzi's state.
******************************************************************* */
int EcoSIM::AdvanceSingleColumn(double dt, int col)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 1;
/*****************************************************************
   ADVANCE CALL GOES HERE
  *******************************************************************/

 bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_sizes_, num_iterations, col);

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyEcoSIMStateToAmanzi(col,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
}

} // namespace EcoSIM
} // namespace Amanzi
