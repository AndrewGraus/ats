/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus


  --------------------------------------------------------------------------*/

#include <algorithm>
#include <set>
#include <string>

// TPLs
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "Mesh.hh"

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
    domain_ = plist_->get<std::string>("domain name", "domain");
    domain_surf_ = Keys::readDomainHint(*plist_, domain_, "subsurface", "surface");

    // obtain key of fields
    // What fields will we need to pass to EcoSIM, presumably fields relating to
    // transport, flow, and energy
    // What do we need for EcoSIM? Based on the variables doc it is:
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
    //
    // transport
    tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");

    //Flow
    poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_liquid_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    saturation_gas_key_ = Keys::readKey(*plist_,domain_,"saturation gas", "saturation_gas");
    saturation_ice_key_ = Keys::readKey(*plist_,domain_,"saturation ice", "saturation_ice");
    //elev_key_ = Keys::readKey(*plist_, domain_, "elevation", "elevation");
    water_content_key_ = Keys::readKey(*plist_,domain_,"water content","water_content");
    rel_perm_key_ = Keys::readKey(*plist_,domain_,"relative permeabiilty","relative_permeability");

    //densities
    //If we need bulk density do we need volume fractions of each quantity?
    //This can be computed from the saturations and porosity (I think) via:
    // f_rock = (1 - porosity)
    // f_liq = S_liq * porosity
    // f_gas = S_gas * porosity
    // f_ice = S_ice * porosity

    liquid_den_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
    ice_den_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    gas_den_key_ = Keys::readKey(*plist_, domain_,"porosity", "porosity");
    rock_den_key_ = Keys::readKey(*plist_, domain_, "density rock", "density_rock");

    //energy
    T_key_ = Keys::readKey(*plist_, domain_, "temperature", "temperature");
    conductivity_key_ = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");

    //Other
    cv_key_ = Keys::readKey(*plist_, domain_, "cell volume", "cell_volume");
    min_vol_frac_key_ = Keys::readKey(*plist_, domain_, "mineral volume fractions", "mineral_volume_fractions");
    ecosim_aux_data_key_ = Keys::readKey(*plist_, domain_, "ecosim aux data", "ecosim_aux_data");

    // parameters
    // initial timestep
    dt_ = plist_->get<double>("initial time step", 1.);
    //Heat capacity looks like the default units are molar heat capacity
    c_m_ = plist_->get<double>("heat capacity [J mol^-1 K^-1]");

    //They also sometimes use a version of heat capacity that is just this
    //quantity times 1e-6:
    //ka_ = 1.e-6 * plist_.get<double>("heat capacity [J kg^-1 K^-1]");
    //Unclear what we need

    // Here is where Alquimia initializes the Chemistry engine which handles differences between
    // CrunchFlow and PFlotran and eventually runs either code to advance the Chemistry
    // We will probably need something like this eventually if we add in other BGC codes but it
    // is very complex and we can almost certainly do something simpler to start out with to just get
    // EcoSIM running.

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

    //comp_names_.clear();
    //bgc_engine_->GetPrimarySpeciesNames(comp_names_);

    //number_aqueous_components_ = comp_names_.size();
    //number_free_ion_ = number_aqueous_components_;
    //number_total_sorbed_ = number_aqueous_components_;

  }

/* *******************************************************************
* Destroy ansilary data structures.
******************************************************************* */
EcoSIM::~EcoSIM()
  {
  if (bgc_initialized_)
    bgc_engine_->FreeState(bgc_props_, bgc_state_, bgc_aux_data_);
  }


// now the PK setup
void EcoSIM::Setup() {
  std::cout << "beginning Ecosim setup\n";
  //PK_Physical_Default::Setup();

  /*This is for setting up the Auxiliary Output data which I'm not sure we need
  chem_engine_->GetAuxiliaryOutputNames(aux_names_, aux_subfield_names_);
  for (size_t i = 0; i < aux_names_.size(); ++i) {
    aux_names_[i] = Keys::getKey(domain_, aux_names_[i]);

    if (!S_->HasRecord(aux_names_[i])) {
      S_->Require<CompositeVector, CompositeVectorSpace>(aux_names_[i], tag_next_, passwd_, aux_subfield_names_[i])
        .SetMesh(mesh_)->SetGhosted(false)
        ->SetComponent("cell", AmanziMesh::CELL, aux_subfield_names_[i].size());
    }
  }
  */

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
  std::cout << "\nEnd setup\n";
}

// -- Initialize owned (dependent) variables.
void EcoSIM::Initialize() {
  std::cout << "\nBegin Initialize\n";
  //PK_Physical_Default::Initialize();

  //Now we have to initalize the variables (i.e. give them initial values)
  //In our PK it will only be done for variables owned by the PK so the aux_names
  //Keeping an example of how it's done generically here:
  //S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_).PutScalar(0.);
  //S_->GetRecordW("co2_decomposition", tag_next_, name_).set_initialized();

  //In alquimia they initialize the axuiliary data via a function called InitializeCVField
  //which can be found in Amanzi/src/PKs/PK_Physical.cc:

  // initialize fields as soon as possible
  /*for (size_t i = 0; i < aux_names_.size(); ++i) {
    InitializeCVField(S_, *vo_, aux_names_[i], tag_next_, passwd_, 0.0);
  }*/

  //Now we call the engine's init state function which allocates the data
  std::cout << "\ninitializing BGC engine with: " << ncells_per_col_ << "\n";

  bgc_engine_->InitState(bgc_props_, bgc_state_, bgc_aux_data_, ncells_per_col_);
  std::cout << "\engine initialized\n";
  //This function calls four separate functions from the interface:
  // AllocateAlquimiaProperties - Allocates the properties which are things
  // chemistry doesn't change
  // AllocateAlquimiaState - Allocates properties of things chemistry CAN change
  // AllocateAqluimiaAuxiliaryData - Allocates variables Alquimia needs to carry
  // over between runs but ATS doesn't need
  // AllocateAlquimiaAuxiliaryOutputData - Allocates variables that ATS will eventually
  // output (probably don't need for now)
  int ierr = 0;

  // Ensure dependencies are fille
  std::cout << "\nfilling dependencies\n";
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(gas_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(conductivity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);
  std::cout << "\ndependencies filled\n";

  int num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //This is the main set up code in alquimia it loops over times and chemical conditions
  //I don't know that we need the two initial loops. I'm just including them because we might
  std::cout << "\ninitializing column loop\n";
  for (int col=0; col!=num_cols_; ++col) {
    //FieldToColumn_(col, temp, col_temp.ptr());
    //ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

    //We're going to need to write an InitializeSingleColumn code
    //ierr = InitializeSingleCell(cell, condition);
    std::cout << "\ninitializing column "<< col <<" \n";
    ierr = InitializeSingleColumn(col);
    //In Alquimia this function simply calls CopyToAlquimia, then
    //Calls the chemistry engine and enforces condition, then copies
    //From Alquimia to Amanzi
    //
    //The copy to alquimia function takes the cell index, because it
    //is assigning things cell by cell in state. For colums this will
    //be a bit trickier I THINK we could use the FieldToColumn_ function
    //To do this, but I think I will actually need to figure out a test
    //for this before I actually code it up
  }
  std::cout << "\nfinishing initialize\n";
  // verbose message
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK was successful, T="
        << S_->get_time() << vo_->reset() << std::endl << std::endl;
  }
}

void EcoSIM::CommitStep(double t_old, double t_new, const Tag& tag) {
  std::cout << "\nRunning commit\n";

  // I don't know that we will have much to do here. In SimpleBGC they just copy
  // Data to the pfts, which we won't be doing. In Alquimia they just save the time
  // As below.

  saved_time_ = t_new;

  std::cout << "\nEnd commit\n";
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
  // Fix this from DEFAULT, see amanzi/amanzi#646 --etc
  // Update all dependencies again
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_liquid_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_gas_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_ice_key_, Tags::DEFAULT).Update(*S_, name_);
  //S_->GetEvaluator(elev_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(water_content_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rel_perm_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(liquid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(gas_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rock_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(conductivity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);

  //---------------------------------------------------------------------------
  //BGCSimple Advance
  //---------------------------------------------------------------------------
  // Copy the PFT from old to new, in case we failed the previous attempt at
  // this timestep.  This is hackery to get around the fact that PFTs are not
  // (but should be) in state.
  AmanziMesh::Entity_ID num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // grab the required fields
  /*
  Epetra_MultiVector& sc_pools = *S_->GetW<CompositeVector>(key_, tag_next_, name_)
      .ViewComponent("cell",false);
  Epetra_MultiVector& co2_decomp = *S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_)
      .ViewComponent("cell",false);
  Epetra_MultiVector& trans = *S_->GetW<CompositeVector>(trans_key_, tag_next_, name_)
      .ViewComponent("cell",false);*/
  //Do I need to update everything here?

  S_->GetEvaluator("porosity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& porosity = *(*S_->Get<CompositeVector>("porosity", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_saturation = *(*S_->Get<CompositeVector>("saturation_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_gas", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& gas_saturation = *(*S_->Get<CompositeVector>("saturation_gas", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("saturation_ice", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& ice_saturation = *(*S_->Get<CompositeVector>("saturation_ice", tag_next_)
          .ViewComponent("cell",false))(0);

  //S_->GetEvaluator("elevation", tag_next_).Update(*S_, name_);
  //const Epetra_MultiVector& elevation = *S_->Get<CompositeVector>("elevation", tag_next_)
  //        .ViewComponent("cell",false);

  S_->GetEvaluator("water_content", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& water_content = *(*S_->Get<CompositeVector>("water_content", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("relatiive_permeabiilty", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& relative_permeability = *(*S_->Get<CompositeVector>("relative_permeabiilty", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("mass_density_liquid", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& liquid_density = *(*S_->Get<CompositeVector>("mass_density_liquid", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("mass_density_gas", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& gas_density = *(*S_->Get<CompositeVector>("mass_density_gas", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("mass_density_ice", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& ice_density = *(*S_->Get<CompositeVector>("mass_density_ice", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("density_rock", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& rock_density = *(*S_->Get<CompositeVector>("density_rock", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& temp = *(*S_->Get<CompositeVector>("temperature", tag_next_)
      .ViewComponent("cell",false))(0);

  S_->GetEvaluator("thermal_conductivity", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& conductivity = *(*S_->Get<CompositeVector>("thermal_conductivity", tag_next_)
          .ViewComponent("cell",false))(0);

  S_->GetEvaluator("cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& cell_volume = *(*S_->Get<CompositeVector>("cell_volume", tag_next_)
          .ViewComponent("cell",false))(0);

  //S_->GetEvaluator("pressure", tag_next_).Update(*S_, name_);
  //const Epetra_MultiVector& pres = *S_->Get<CompositeVector>("pressure", tag_next_)
  //    .ViewComponent("cell",false);

  // note that this is used as the column area, which is maybe not always
  // right.  Likely correct for soil carbon calculations and incorrect for
  // surface vegetation calculations (where the subsurface's face area is more
  // correct?)
  S_->GetEvaluator("surface-cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& scv = *(*S_->Get<CompositeVector>("surface-cell_volume", tag_next_)
      .ViewComponent("cell", false))(0);

  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col=0; col!=num_cols_; ++col) {

    auto& col_iter = mesh_->cells_of_column(col);
    ncells_per_col_ = col_iter.size();

    //Copy to EcoSIM structures

    // call the model
    // This will be where we call the main function which advances ecosim
    //BGCAdvance(S_->get_time(tag_current_), dt, scv[0][col], cryoturbation_coef_, met,
    //           *temp_c, *pres_c, *depth_c, *dz_c,
    //           pfts_[col], soil_carbon_pools_[col],
    //           co2_decomp_c, trans_c, sw_c);

    AdvanceSingleColumn(dt, col);

    //Copy back to Amanzi

  } // end loop over columns

  // mark primaries as changed

  // Compute the next time step.
  // * will we need to do this? *
  //ComputeNextTimeStep();

  //return failed;

  std::cout << "\nEnd Advance\n";
}

//---------------------------------------------------------------------------
//Here are the BGCSimple helper functions
//---------------------------------------------------------------------------

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
       Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  std::cout << "\ncol: "<< col <<"\n";
  auto& col_iter = mesh_->cells_of_column(col);
  std::cout << "\ncol_iter: "<< col_iter.size() <<"\n";

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    std::size_t vec_index = col_iter[i];
    std::cout << "for i: " << i << "vec_index: " << vec_index << "\n";
    std::cout << "vec[" << vec_index << "]: " << vec[vec_index] << "\n";
  }


  for (std::size_t i=0; i!=col_iter.size(); ++i) {

    if (i >= col_vec->Length()) {
      std::cout << "Error: index " << i << " is out of bounds for col_vec\n";
    }

    std::size_t vec_index = col_iter[i];

    if (vec_index >= vec.MyLength()) {
      std::cout << "Error: index " << vec_index << " is out of bounds for vec\n";
    }

    std::cout << "col_vec[" << i << "]: " << (*col_vec)[i] << "\n";
    std::cout << "vec[" << vec_index << "]: " << vec[vec_index] << "\n";

    //(*col_vec)[i] = vec[vec_index];
    (*col_vec)[i] = vec[vec_index];
  }
}

// helper function for pushing field to column
/*void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_MultiVector& vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    col_vec[i] = vec[col_iter[i]];
  }
}*/

// I think I need a function for pushing from the column back to the field
// with any luck it's just the reverse of the above similar to how it's done
// cell by cell in alquimia

void EcoSIM::ColumnToField_(AmanziMesh::Entity_ID col, Epetra_Vector& vec,
                               Teuchos::Ptr<Epetra_SerialDenseVector> col_vec)
{
  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    vec[col_iter[i]] = (*col_vec)[i];
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

//---------------------------------------------------------------------------
//Alquimia Helper functions
//---------------------------------------------------------------------------
/*void EcoSIM::CopyToEcoSIM(int col,
        BGCProperties& props,
        BGCState& state,
        BGCAuxiliaryData& aux_data)
{
  CopyToEcoSIM(col, props, state, aux_data);
}*/

void EcoSIM::CopyToEcoSIM(int col,
                                 BGCProperties& props,
                                 BGCState& state,
                                 BGCAuxiliaryData& aux_data,
                               const Tag& water_tag)
{
  //Fill state with ATS variables that are going to be changed by EcoSIM
  //NEED TO DECIDE WHICH PROPERTIES GO WHERE
  std::cout << "\nviewing components\n";
  //Might need to switch
  //const Epetra_MultiVector& temp also set ViewComponent to false
  const Epetra_Vector& porosity = *(*S_->Get<CompositeVector>(poro_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& tcc = *(*S_->Get<CompositeVector>(tcc_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_saturation = *(*S_->Get<CompositeVector>(saturation_liquid_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& gas_saturation = *(*S_->Get<CompositeVector>(saturation_gas_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& ice_saturation = *(*S_->Get<CompositeVector>(saturation_ice_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& elevation = *S_->Get<CompositeVector>(elev_key_, water_tag).ViewComponent("cell", false);
  const Epetra_Vector& water_content = *(*S_->Get<CompositeVector>(water_content_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& relative_permeability = *(*S_->Get<CompositeVector>(rel_perm_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& liquid_density = *(*S_->Get<CompositeVector>(liquid_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& ice_density = *(*S_->Get<CompositeVector>(ice_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& gas_density = *(*S_->Get<CompositeVector>(gas_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& rock_density = *(*S_->Get<CompositeVector>(rock_den_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& temp = *(*S_->Get<CompositeVector>(T_key_, water_tag).ViewComponent("cell", false))(0);
  //const Epetra_Vector& temp = *S_->Get<CompositeVector>("temperature", tag_next_).ViewComponent("cell", false);

  const Epetra_Vector& conductivity = *(*S_->Get<CompositeVector>(conductivity_key_, water_tag).ViewComponent("cell", false))(0);
  const Epetra_Vector& cell_volume = *(*S_->Get<CompositeVector>(cv_key_, water_tag).ViewComponent("cell", false))(0);

  //Define the column vectors to hold the data
  std::cout << "\ncreating column vectors with size: "<< ncells_per_col_ << "\n";
  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  std::cout << "\ncreated first column\n";
  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_elev = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_f_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  //Here is where we should do the various field-to-column calls to then pass along
  //to the data structures that will pass the data to EcoSIM
  //Format is:
  //FieldToColumn_(column index, dataset to copy from, vector to put the data in)

  //FieldToColumn_(col,tcc,col_tcc.ptr());
  std::cout << "\nCopying Amanzi field to column vector\n";
  FieldToColumn_(col,porosity,col_poro.ptr());
  std::cout << "\npushed first column\n";
  FieldToColumn_(col,liquid_saturation,col_l_sat.ptr());
  FieldToColumn_(col,gas_saturation,col_g_sat.ptr());
  FieldToColumn_(col,ice_saturation,col_i_sat.ptr());
  //FieldToColumn_(col,elevation,col_elev.ptr());
  FieldToColumn_(col,water_content,col_wc.ptr());
  FieldToColumn_(col,relative_permeability,col_rel_perm.ptr());
  FieldToColumn_(col,liquid_density,col_f_dens.ptr());
  FieldToColumn_(col,ice_density,col_i_dens.ptr());
  FieldToColumn_(col,gas_density,col_g_dens.ptr());
  FieldToColumn_(col,rock_density,col_r_dens.ptr());
  FieldToColumn_(col,temp, col_temp.ptr());
  FieldToColumn_(col,conductivity,col_cond.ptr());
  FieldToColumn_(col,cell_volume,col_vol.ptr());

  // I think I need to loop over the column data and save it to the data
  // structures. Eventually I could probably rewrite FieldToColumn_ to do this
  // automatically, but I just want to test this for now

  std::cout << "\nlooping over cells and copying to EcoSIM data structure\n";
  for (int i=0; i < ncells_per_col_; ++i) {
    std::cout << "\nlooping through cell " << i << "\n";
    state.fluid_density.data[i] = (*col_f_dens)[i];
    state.gas_density.data[i] = (*col_g_dens)[i];
    state.ice_density.data[i] = (*col_i_dens)[i];
    state.porosity.data[i] = (*col_poro)[i];
    state.water_content.data[i] = (*col_wc)[i];
    state.temperature.data[i] = (*col_temp)[i];

    props.liquid_saturation.data[i] = (*col_l_sat)[i];
    props.gas_saturation.data[i] = (*col_g_sat)[i];
    props.ice_saturation.data[i] = (*col_i_sat)[i];
    //props.elevation.data[i] = (*col_elev)[i];
    props.relative_permeability.data[i] = (*col_rel_perm)[i];
    props.conductivity.data[i] = (*col_cond)[i];
    props.volume.data[i] = (*col_vol)[i];

  }
  std::cout << "\nFinished loop\n";
  //mat_props.volume = mesh_->cell_volume(cell;z
  //mat_props.saturation = water_saturation[0][cell];

  num_components_ = tcc.NumVectors();

  //This probably isn't going to work. I think I either need to think
  //of a way to do this
  //Possible ideas:
  // 1) creat a data column per component (how to do that without hard coding?)
  // 2) have some sort of array of shape num_componentsXcolumnsize and loop over
  //    the components
  // For #2 I think I just need to change the serieal dense vector call to
  // a different data type (are these always 1d?) What is the 2d version?
  //for (int i = 0; i < num_components; i++) {
  //  bgc_state.total_mobile.data[i] = (*col_tcc)[i];
  //}

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
  //I probably need to copy the columns cell by cell in a loop
  //Can I do this in the field to column function?

  //auto& tcc = S_->GetPtrW<CompositeVector>(tcc_key_, water_tag, passwd_).ViewComponent("cell");

  //There are various ways to access data unclear which I need
  //Attempt 1 fails
  //Epetra_MultiVector&  porosity = *S_->GetPtrW<CompositeVector>(poro_key_, water_tag, passwd_)->ViewComponent("cell");

  //Attempt 2
  //porosity = Teuchos::rcp(new Epetra_MultiVector(*S_->Get<CompositeVector>("porosity").ViewComponent("cell", true)));

  //Attempt 3
  //porosity = S_->GetPtrW<CompositeVector>(poro_key_, water_tag, passwd_)->ViewComponent("cell", false));

  //Attempt 4
  //Teuchos::RCP<Epetra_MultiVector> porosity = S_->GetW<CompositeVector>(poro_key_, water_tag, passwd_).ViewComponent("cell", true);

  //Attempt 5 (ELM; this should work)
  auto& porosity = *(*S_->GetW<CompositeVector>(poro_key_, Amanzi::Tags::NEXT, poro_key_).ViewComponent("cell",false))(0);

  auto& liquid_saturation = *(*S_->GetW<CompositeVector>(saturation_liquid_key_, Amanzi::Tags::NEXT, saturation_liquid_key_).ViewComponent("cell",false))(0);
  auto& gas_saturation = *(*S_->GetW<CompositeVector>(saturation_gas_key_, Amanzi::Tags::NEXT, saturation_gas_key_).ViewComponent("cell",false))(0);
  auto& ice_saturation = *(*S_->GetW<CompositeVector>(saturation_ice_key_, Amanzi::Tags::NEXT, saturation_ice_key_).ViewComponent("cell",false))(0);
  //auto& elevation = S_->GetPtrW<CompositeVector>(elev_key_, Amanzi::Tags::NEXT, passwd_).ViewComponent("cell");
  auto& water_content = *(*S_->GetW<CompositeVector>(water_content_key_, Amanzi::Tags::NEXT, water_content_key_).ViewComponent("cell",false))(0);
  auto& relative_permeability = *(*S_->GetW<CompositeVector>(rel_perm_key_, Amanzi::Tags::NEXT, rel_perm_key_).ViewComponent("cell",false))(0);
  auto& liquid_density = *(*S_->GetW<CompositeVector>(liquid_den_key_, Amanzi::Tags::NEXT, liquid_den_key_).ViewComponent("cell",false))(0);
  auto& ice_density = *(*S_->GetW<CompositeVector>(ice_den_key_, Amanzi::Tags::NEXT, ice_den_key_).ViewComponent("cell",false))(0);
  auto& gas_density = *(*S_->GetW<CompositeVector>(gas_den_key_, Amanzi::Tags::NEXT, gas_den_key_).ViewComponent("cell",false))(0);
  auto& rock_density = *(*S_->GetW<CompositeVector>(rock_den_key_, Amanzi::Tags::NEXT, rock_den_key_).ViewComponent("cell",false))(0);
  //auto& temp = *S_->GetW<CompositeVector>(T_key_, Amanzi::Tags::NEXT, name()).ViewComponent("cell",false);
  auto& conductivity = *(*S_->GetW<CompositeVector>(conductivity_key_, Amanzi::Tags::NEXT, conductivity_key_).ViewComponent("cell",false))(0);
  auto& cell_volume = *(*S_->GetW<CompositeVector>(cv_key_, Amanzi::Tags::NEXT, cv_key_).ViewComponent("cell",false))(0);

  //I think I need to redefine this here?
  auto col_tcc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_poro = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_l_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_sat = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_elev = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_wc = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_rel_perm = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_f_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_i_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_g_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_r_dens = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  //auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_cond = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_vol = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  /*Can't save cell by cell doesn't seem to work like this*/

  for (int i=0; i < ncells_per_col_; ++i) {
    (*col_f_dens)[i] = state.fluid_density.data[i];
    (*col_g_dens)[i] = state.gas_density.data[i];
    (*col_i_dens)[i] = state.ice_density.data[i];
    (*col_poro)[i] = state.porosity.data[i];
    (*col_wc)[i] = state.water_content.data[i];
    //(*col_temp)[i] = state.temperature.data[i];

    (*col_l_sat)[i] = props.liquid_saturation.data[i];
    (*col_g_sat)[i] = props.gas_saturation.data[i];
    (*col_i_sat)[i] = props.ice_saturation.data[i];
    //(*col_elev)[i] = props.elevation.data[i];
    (*col_rel_perm)[i] = props.relative_permeability.data[i];
    (*col_cond)[i] = props.conductivity.data[i];
    (*col_vol)[i] = props.volume.data[i];
  }

  /*for (int i = 0; i < num_components; i++) {
    bgc_state.total_mobile.data[i] = (*col_tcc)[i];
  }*/

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
  ColumnToField_(col,gas_saturation,col_g_sat.ptr());
  ColumnToField_(col,ice_saturation,col_i_sat.ptr());
  //ColumnToField_(col,elevation,col_elev.ptr());
  ColumnToField_(col,water_content,col_wc.ptr());
  ColumnToField_(col,relative_permeability,col_rel_perm.ptr());
  ColumnToField_(col,liquid_density,col_f_dens.ptr());
  ColumnToField_(col,ice_density,col_i_dens.ptr());
  ColumnToField_(col,gas_density,col_g_dens.ptr());
  ColumnToField_(col,rock_density,col_r_dens.ptr());
  //ColumnToField_(col,temp, col_temp.ptr());
  ColumnToField_(col,conductivity,col_cond.ptr());
  ColumnToField_(col,cell_volume,col_vol.ptr());

}

/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int EcoSIM::InitializeSingleColumn(int col)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  std::cout << "\nCopying to EcoSIM\n";
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);
  std::cout << "\nFinished Copying to EcoSIM\n";

  //bgc_engine_->EnforceCondition(condition, current_time_, bgc_props_,
  //        bgc_state_, bgc_aux_data_);
  std::cout << "\nCopying back to Amanzi\n";
  CopyEcoSIMStateToAmanzi(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);
  std::cout << "\nFinished copying to Amanzi\n";

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
* It returns the number of iterations taken to obtain the advanced solution,
* or -1 if an error occurred.
******************************************************************* */
int EcoSIM::AdvanceSingleColumn(double dt, int col)
{
  // Copy the state and property information from Amanzi's state within
  // this cell to Alquimia.
  //
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  int num_iterations = 0;

  //Think a bit about what to do with this
  /*****************************************************************
   ADVANCE CALL GOES HERE
  ******************************************************************

  if (ecosim_mat_props_.saturation > saturation_tolerance_) {
    bool success = bgc_engine_->Advance(dt, bgc_props_, bgc_state_,
                                         bgc_aux_data_, num_iterations);
    if (not success) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "no convergence in cell: " << mesh_->cell_map(false).GID(cell) << std::endl;
      }
      return -1;
    }
  }
  */

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyEcoSIMStateToAmanzi(col,
                            bgc_props_, bgc_state_, bgc_aux_data_, Tags::DEFAULT);

  return num_iterations;
}

} // namespace EcoSIM
} // namespace Amanzi
