/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus


  --------------------------------------------------------------------------*/

#include "Ecosim_ATS_interface.hh"
#include "pk_helpers.hh"

namespace Amanzi {
namespace EcoSIM {

EcoSIM::EcoSIM(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution):
  PK_Physical_Default(pk_tree, global_list, S, solution),
  PK(pk_tree, global_list, S, solution)
  {
    domain_ = plist_->get<std::string>("domain name", "domain");

    // obtain key of fields
    // What fields will we need to pass to EcoSIM, presumably fields relating to
    // transport, flow, and energy
    // What do we need for EcoSIM? Based on the variables doc it is:
    // grid position (X,Y,Z) - can we just pass Z and assume X and Y are 0?
    // soil texture - not sure where this is
    // bulk density - not sure
    // Elevation
    // Aspect in geometric format
    // soil texture (sand, clay, silt)
    // water table depth

    // transport
    tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");

    //Flow
    poro_key_ = Keys::readKey(*plist_, domain_, "porosity", "porosity");
    saturation_key_ = Keys::readKey(*plist_, domain_, "saturation liquid", "saturation_liquid");
    fluid_den_key_ = Keys::readKey(*plist_, domain_, "mass density liquid", "mass_density_liquid");
    ice_den_key_ = Keys::readKey(plist_, domain, "ice mass density", "mass_density_ice");
    mass_den_key_ = Keys::readKey(plist_, domain, "mass density", mass_key);
    rhos_key_ = Keys::readKey(plist_, domain_name, "density rock", "density_rock");

    //energy
    T_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
    conductivity_key_ = = Keys::readKey(*plist_, domain_, "thermal conductivity", "thermal_conductivity");

    //Other
    cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
    min_vol_frac_key_ = Keys::readKey(*plist_, domain_, "mineral volume fractions", "mineral_volume_fractions");
    ecosim_aux_data_key_ = Keys::readKey(*plist_, domain_, "ecosim aux data", "ecosim_aux_data");

    // parameters
    // initial timestep
    dt_ = plist_->get<double>("initial time step", 1.);
    //Heat capacity looks like the default units are molar heat capacity
    c_m = plist_.get<double>("heat capacity [J mol^-1 K^-1]");
    //They also sometimes use a version of heat capacity that is just this
    //quantity times 1e-6:
    //ka_ = 1.e-6 * plist_.get<double>("heat capacity [J kg^-1 K^-1]");
    //Unclear what we need

    // Here is where Alquimia initializes the Chemistry engine which handles differences between
    // CrunchFlow and PFlotran and eventually runs either code to advance the Chemistry
    // We will probably need something like this eventually if we add in other BGC codes but it
    // is very complex and we can almost certainly do something simpler to start out with to just get
    // EcoSIM running.
    //
    // For now I'm changing the syntax from chem_engine to bgc_engine but I'll probably end up
    // either cutting this out or commenting it out.

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
    bgc_engine_ = Teuchos::rcp(new EcoSIM::BGCEngine(engine_name, engine_inputfile));
    //When creating the engine there are a few functions in ChemistryEngine that
    //are called AllocateAlquimiaEngineStatus, CreateAlquimiaInterface, and chem_.Setup
    //these are functions in the Alquimia standalone code
    //
    // AllocateAlquimiaEngineStatus - is in alquimia_memory.c does some simple memory
    // Calculation I think.
    //
    // CreateAlquimiaInterface - in alquimia_interface.c, this makes the main
    // decision between Pflotran and CrunchFlow, then allocates to the interface
    // called chem_ where to find the processes the interface will need including
    // Setup, Shutdown, ProcessCondition, ReactionStepOperatorSplit, GetAuxiliaryOutput
    // GetProblemMetaData
    //
    // Basically it goes Alquimia PK -> ChemistryEngine -> Alquimia Interface ->
    // CrunchFlow/PFloTran. We don't have to stick to this and probably just simplify to:
    // EcoSIM_PK -> BGCEngine -> EcoSIM driver

    // grab the component names
    // This is a small function in ChemistryEngine that simply looks at the metadata
    // finds the number of primary species and fills a vector with the names of those
    //species
    comp_names_.clear();
    bgc_engine_->GetPrimarySpeciesNames(comp_names_);

    number_aqueous_components_ = comp_names_.size();
    number_free_ion_ = number_aqueous_components_;
    number_total_sorbed_ = number_aqueous_components_;

  }

/* *******************************************************************
* Destroy ansilary data structures.
******************************************************************* */
//Going to need this eventually to clear the strutures
//the various data structures were named for alquimia i.e.
//alq_mat_props_, alq_state, ect. changed to just bgc
EcoSIM::~EcoSIM()
  {
  if (bgc_initialized_)
    bgc_engine_->FreeState(bgc_props_, bgc_state_, bgc_aux_data_);
  }


// now the PK setup
void EcoSIM::Setup() {
  std::cout << "beginning Ecosim setup\n";
  PK_Physical_Default::Setup();
  //I actually don't think we have much to do here. In BGC_simple they set up
  //The arrays for the PFTs and the carbon pools but we're not going to setup
  //variables in this way. Other than that it sets the required variables. It
  //seems like these are things that are "owned" by this PK so the auxiliary data
  //Leaving in the alquimia code for this. It seems to set the auxiliary data up
  //in several different ways. Still need to understand what is going on here.
  //
  //To further confuse things the chemistry engine function called Setup occurs
  //When the engine is created in the constructor NOT called in the PK setup

  // Set up auxiliary chemistry data using the ChemistryEngine.

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

  //This is for the Auxiliary Data which we will need
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
    int num_aux_data = chem_engine_->Sizes().num_aux_integers + chem_engine_->Sizes().num_aux_doubles;
    S_->Require<CompositeVector, CompositeVectorSpace>(alquimia_aux_data_key_, tag_next_, passwd_)
      .SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, num_aux_data);

    S_->GetRecordW(alquimia_aux_data_key_, tag_next_, passwd_).set_io_vis(false);
  }

  std::cout << "\nEnd setup\n";
}

// -- Initialize owned (dependent) variables.
void EcoSIM::Initialize() {
  std::cout << "\nBegin Initialize\n";
  //PK_Physical_Default::Initialize(S);
  PK_Physical_Default::Initialize();

  //Now we have to initalize the variables (i.e. give them initial values)
  //In our PK it will only be done for variables owned by the PK so the aux_names
  //Keeping an example of how it's done generically here:
  //S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_).PutScalar(0.);
  //S_->GetRecordW("co2_decomposition", tag_next_, name_).set_initialized();

  //In alquimia they initialize the axuiliary data via a function called InitializeCVField
  //which can be found in Amanzi/src/PKs/PK_Physical.cc:

  // initialize fields as soon as possible
  for (size_t i = 0; i < aux_names_.size(); ++i) {
    InitializeCVField(S_, *vo_, aux_names_[i], tag_next_, passwd_, 0.0);
  }

  //Now we call the engine's init state function which allocates the data
  bgc_engine_->InitState(bgc_props_, bgc_state_, bgc_aux_data_);
  //This function calls four separate functions from the interface:
  // AllocateAlquimiaProperties - Allocates the properties which are things
  // chemistry doesn't change
  // AllocateAlquimiaState - Allocates properties of things chemistry CAN change
  // AllocateAqluimiaAuxiliaryData - Allocates variables Alquimia needs to carry
  // over between runs but ATS doesn't need
  // AllocateAlquimiaAuxiliaryOutputData - Allocates variables that ATS will eventually
  // output (probably don't need for now)

  // Ensure dependencies are filled
  // I think this is everything from ATS that needs to be updated
  S_->GetEvaluator(tcc_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(poro_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(fluid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(mass_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rhos_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(T_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(conductivity_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(cv_key_, Tags::DEFAULT).Update(*S_, name_);

  // init root carbon
  // I think we can basically make use of the FieldToColumn_ function to take
  // all of the properties we need (water, heat, ect) and pass them to each
  // column
  auto col_temp = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_depth = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));
  auto col_dz = Teuchos::rcp(new Epetra_SerialDenseVector(ncells_per_col_));

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_Vector& temp = *(*S_->Get<CompositeVector>("temperature", tag_next_)
				.ViewComponent("cell",false))(0);

  int num_cols_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  //This is the main set up code in alquimia it loops over times and chemical conditions
  //I don't know that we need the two initial loops. I'm just including them because we might

  if (fabs(initial_conditions_time_ - S_->get_time()) < 1e-8 * (1.0 + fabs(S_->get_time()))) {
    for (auto it = chem_initial_conditions_.begin(); it != chem_initial_conditions_.end(); ++it) {
      std::string region = it->first;
      std::string condition = it->second;

      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "enforcing geochemical condition \"" << condition
                   << "\" in region \"" << region << "\"\n";
            }
            for (int col=0; col!=num_cols_; ++col) {
              FieldToColumn_(col, temp, col_temp.ptr());
              ColDepthDz_(col, col_depth.ptr(), col_dz.ptr());

              //We're going to need to write an InitializeSingleColumn code
              //ierr = InitializeSingleCell(cell, condition);

              ierr = InitializeSingleColumn(col, condition);
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
  S_->GetEvaluator(fluid_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(saturation_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(ice_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(mass_den_key_, Tags::DEFAULT).Update(*S_, name_);
  S_->GetEvaluator(rhos_key_, Tags::DEFAULT).Update(*S_, name_);
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
  Epetra_MultiVector& sc_pools = *S_->GetW<CompositeVector>(key_, tag_next_, name_)
      .ViewComponent("cell",false);
  Epetra_MultiVector& co2_decomp = *S_->GetW<CompositeVector>("co2_decomposition", tag_next_, name_)
      .ViewComponent("cell",false);
  Epetra_MultiVector& trans = *S_->GetW<CompositeVector>(trans_key_, tag_next_, name_)
      .ViewComponent("cell",false);

  S_->GetEvaluator("temperature", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& temp = *S_->Get<CompositeVector>("temperature", tag_next_)
      .ViewComponent("cell",false);

  S_->GetEvaluator("pressure", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& pres = *S_->Get<CompositeVector>("pressure", tag_next_)
      .ViewComponent("cell",false);

  // note that this is used as the column area, which is maybe not always
  // right.  Likely correct for soil carbon calculations and incorrect for
  // surface vegetation calculations (where the subsurface's face area is more
  // correct?)
  S_->GetEvaluator("surface-cell_volume", tag_next_).Update(*S_, name_);
  const Epetra_MultiVector& scv = *S_->Get<CompositeVector>("surface-cell_volume", tag_next_)
      .ViewComponent("cell", false);


  // loop over columns and apply the model
  for (AmanziMesh::Entity_ID col=0; col!=num_cols_; ++col) {
    // update the various soil arrays
    FieldToColumn_(col, *temp(0), temp_c.ptr());
    FieldToColumn_(col, *pres(0), pres_c.ptr());
    ColDepthDz_(col, depth_c.ptr(), dz_c.ptr());

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
  changedEvaluatorPrimary(trans_key_, tag_next_, *S_);
  changedEvaluatorPrimary(shaded_sw_key_, tag_next_, *S_);
  changedEvaluatorPrimary(total_lai_key_, tag_next_, *S_);
  return false;

  // Compute the next time step.
  // * will we need to do this? *
  ComputeNextTimeStep();

  return failed;

  std::cout << "\nEnd Advance\n";
}

//---------------------------------------------------------------------------
//Here are the BGCSimple helper functions
//---------------------------------------------------------------------------

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
        Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy)
{
  if (col_vec == Teuchos::null) {
    col_vec = Teuchos::ptr(new Epetra_SerialDenseVector(ncells_per_col_));
  }

  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    (*col_vec)[i] = vec[col_iter[i]];
  }
}

// helper function for pushing field to column
void EcoSIM::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                               double* col_vec, int ncol)
{
  auto& col_iter = mesh_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    col_vec[i] = vec[col_iter[i]];
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
void EcoSIM::CopyToEcoSIM(int col,
        BGCProperties& props,
        BGCState& state,
        BGCAuxiliaryData& aux_data)
{
  CopyToEcoSIM(col, props, state, aux_data);
}

void EcoSIM::CopyToEcoSIM(int col,
                                 BGCProperties& props,
                                 BGCState& state,
                                 BGCAuxiliaryData& aux_data)
{
  //NOTE: I Have not touched this yet as it will depend on what we have to
  //Transfer
  //
  //Fill state with ATS variables that are going to be changed by EcoSIM
  const auto& porosity = *S_->Get<CompositeVector>(poro_key_, water_tag).ViewComponent("cell", true);
  const auto& fluid_density = *S_->Get<CompositeVector>(fluid_den_key_, water_tag).ViewComponent("cell", true);
  const auto& water_saturation = *S_->Get<CompositeVector>(saturation_key_, water_tag).ViewComponent("cell", true);
  const auto& tcc = *S_->Get<CompositeVector>(tcc_key_, water_tag).ViewComponent("cell", true);
  const auto& ice_density = *S_->Get<CompositeVector>(ice_den_key_, water_tag).ViewComponent("cell", true);
  const auto& mass_density = *S_->Get<CompositeVector>(mass_den_key_, water_tag).ViewComponent("cell", true);
  const auto& rock_density = *S_->Get<CompositeVector>(rhos_key_, water_tag).ViewComponent("cell", true);
  const auto& Temp = *S_->Get<CompositeVector>(T_key_, water_tag).ViewComponent("cell", true);
  const auto& conductivity = *S_->Get<CompositeVector>(conductivity_key_, water_tag).ViewComponent("cell", true);
  const auto& cell_volume = *S_->Get<CompositeVector>(cv_key_, water_tag).ViewComponent("cell", true);

  state.water_density = fluid_density[0][cell];
  state.porosity = porosity[0][cell];

  //We are probably going to need to do something like this that loops over
  //All transport components to give them over to EcoSIM
  for (int i = 0; i < number_aqueous_components_; i++) {
    state.total_mobile.data[i] = (*aqueous_components)[i][cell];

    if (using_sorption_) {
      const auto& sorbed = *S_->Get<CompositeVector>(total_sorbed_key_, tag_next_).ViewComponent("cell");
      state.total_immobile.data[i] = sorbed[i][cell];
    }
  }

  // minerals
  assert(state.mineral_volume_fraction.size == number_minerals_);
  assert(state.mineral_specific_surface_area.size == number_minerals_);
  assert(mat_props.mineral_rate_cnst.size == number_minerals_);

  if (number_minerals_ > 0) {
    const auto& mineral_vf = *S_->Get<CompositeVector>(min_vol_frac_key_, tag_next_).ViewComponent("cell");
    const auto& mineral_ssa = *S_->Get<CompositeVector>(min_ssa_key_, tag_next_).ViewComponent("cell");
    const auto& mineral_rate = *S_->Get<CompositeVector>(mineral_rate_constant_key_, tag_next_).ViewComponent("cell");
    for (unsigned int i = 0; i < number_minerals_; ++i) {
      state.mineral_volume_fraction.data[i] = mineral_vf[i][cell];
      mat_props.mineral_rate_cnst.data[i] = mineral_rate[i][cell];
      state.mineral_specific_surface_area.data[i] = mineral_ssa[i][cell];
    }
  }

  // ion exchange
  assert(state.cation_exchange_capacity.size == number_ion_exchange_sites_);
  if (number_ion_exchange_sites_ > 0) {
    const auto& ion_exchange = *S_->Get<CompositeVector>(ion_exchange_sites_key_, tag_next_).ViewComponent("cell");
    for (int i = 0; i < number_ion_exchange_sites_; i++) {
      state.cation_exchange_capacity.data[i] = ion_exchange[i][cell];
    }
  }

  // surface complexation
  if (number_sorption_sites_ > 0) {
    const auto& sorption_sites = *S_->Get<CompositeVector>(sorp_sites_key_, tag_next_).ViewComponent("cell");

    assert(number_sorption_sites_ == state.surface_site_density.size);
    for (int i = 0; i < number_sorption_sites_; ++i) {
      // FIXME: Need site density names, too?
      state.surface_site_density.data[i] = sorption_sites[i][cell];
      // TODO(bandre): need to save surface complexation free site conc here!
    }
  }

  // Auxiliary data -- block copy.
  if (S_->HasRecord(alquimia_aux_data_key_, tag_next_)) {
    aux_data_ = S_->GetW<CompositeVector>(alquimia_aux_data_key_, tag_next_, passwd_).ViewComponent("cell");
    int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      aux_data.aux_ints.data[i] = (int)cell_aux_ints[cell];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      aux_data.aux_doubles.data[i] = cell_aux_doubles[cell];
    }
  }

  mat_props.volume = mesh_->cell_volume(cell);
  mat_props.saturation = water_saturation[0][cell];

  // sorption isotherms
  if (using_sorption_isotherms_) {
    const auto& isotherm_kd = *S_->Get<CompositeVector>(isotherm_kd_key_, tag_next_).ViewComponent("cell");
    const auto& isotherm_freundlich_n = *S_->Get<CompositeVector>(isotherm_freundlich_n_key_, tag_next_).ViewComponent("cell");
    const auto& isotherm_langmuir_b = *S_->Get<CompositeVector>(isotherm_langmuir_b_key_, tag_next_).ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      mat_props.isotherm_kd.data[i] = isotherm_kd[i][cell];
      mat_props.freundlich_n.data[i] = isotherm_freundlich_n[i][cell];
      mat_props.langmuir_b.data[i] = isotherm_langmuir_b[i][cell];
    }
  }

  // first order reaction rate cnst
  if (number_aqueous_kinetics_ > 0) {
    const auto& aqueous_kinetics_rate = *S_->Get<CompositeVector>(first_order_decay_constant_key_, tag_next_).ViewComponent("cell");
    for (unsigned int i = 0; i < number_aqueous_kinetics_; ++i) {
      mat_props.aqueous_kinetic_rate_cnst.data[i] = aqueous_kinetics_rate[i][cell];
    }
  }
}

void EcoSIM::CopyEcoSIMStateToAmanzi(
    const int cell,
    const BGCProperties& props,
    const BGCState& state,
    const BGCAuxiliaryData& aux_data)
{
  CopyFromEcoSIM(cell, props, state, aux_data);
}

void EcoSIM::CopyFromEcoSIM(const int cell,
                                   const BGCProperties& props,
                                   const BGCState& state,
                                   const BGCAuxiliaryData& aux_data)
{
  // If the chemistry has modified the porosity and/or density, it needs to
  // be updated here.
  // (this->water_density())[cell] = state.water_density;
  // (this->porosity())[cell] = state.porosity;

  for (int i = 0; i < number_aqueous_components_; ++i) {
    (*aqueous_components)[i][cell] = state.total_mobile.data[i];

    if (using_sorption_) {
      auto& sorbed = *S_->GetW<CompositeVector>(total_sorbed_key_, tag_next_, passwd_).ViewComponent("cell");
      sorbed[i][cell] = state.total_immobile.data[i];
    }
  }

  //Here is where the auxiliary data is filled need to try to change this to columns
  //This may not be trivial
  if (S_->HasRecord(alquimia_aux_data_key_, tag_next_)) {
    aux_data_ = S_->GetW<CompositeVector>(alquimia_aux_data_key_, tag_next_, passwd_).ViewComponent("cell");

    int num_aux_ints = chem_engine_->Sizes().num_aux_integers;
    int num_aux_doubles = chem_engine_->Sizes().num_aux_doubles;

    for (int i = 0; i < num_aux_ints; i++) {
      double* cell_aux_ints = (*aux_data_)[i];
      cell_aux_ints[cell] = (double)aux_data.aux_ints.data[i];
    }
    for (int i = 0; i < num_aux_doubles; i++) {
      double* cell_aux_doubles = (*aux_data_)[i + num_aux_ints];
      cell_aux_doubles[cell] = aux_data.aux_doubles.data[i];
    }
  }


  // Here is where constants are saved as properties (things that won't be
  //changed by EcoSIM)
  if (using_sorption_isotherms_) {
    auto& isotherm_kd = *S_->GetW<CompositeVector>(isotherm_kd_key_, tag_next_, passwd_).ViewComponent("cell");
    auto& isotherm_freundlich_n = *S_->GetW<CompositeVector>(isotherm_freundlich_n_key_, tag_next_, passwd_).ViewComponent("cell");
    auto& isotherm_langmuir_b = *S_->GetW<CompositeVector>(isotherm_langmuir_b_key_, tag_next_, passwd_).ViewComponent("cell");

    for (unsigned int i = 0; i < number_aqueous_components_; ++i) {
      isotherm_kd[i][cell] = mat_props.isotherm_kd.data[i];
      isotherm_freundlich_n[i][cell] = mat_props.freundlich_n.data[i];
      isotherm_langmuir_b[i][cell] = mat_props.langmuir_b.data[i];
    }
  }
}

/* *******************************************************************
* This helper performs initialization on a single cell within Amanzi's state.
* It returns an error code that indicates success (0) or failure (1).
******************************************************************* */
int EcoSIM::InitializeSingleColumn(int col, const std::string& condition)
{
  // NOTE: this should get set not to be hard-coded to Tags::DEFAULT, but
  // should use the same tag as transport.  See #673
  CopyToAlquimia(col, bgc_props_, bgc_state_, bgc_aux_data_);

  bgc_engine_->EnforceCondition(condition, current_time_, bgc_props_,
          bgc_state_, bgc_aux_data_);

  CopyAlquimiaStateToAmanzi(col, bgc_props_, bgc_state_, bgc_aux_data_);


  // ETC: hacking to get consistent solution -- if there is no water
  // (e.g. surface system, we still need to call EnforceCondition() as it also
  // gets aux data set up correctly.  But the concentrations need to be
  // overwritten as 0 to get expected output.  Therefore we manually overwrite
  // this now.  Previously this happened due to a bug in ATS's reactive
  // transport coupler -- happy accidents.
  if (alq_mat_props_.saturation <= saturation_tolerance_)
    for (int i=0; i!=aqueous_components_->NumVectors(); ++i) (*aqueous_components_)[i][cell] = 0.;
  return 0;
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
  CopyToEcoSIM(col, bgc_props_, bgc_state_, bgc_aux_data_);

  int num_iterations = 0;

  //Think a bit about what to do with this
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

  // Move the information back into Amanzi's state, updating the given total concentration vector.
  CopyEcoSIMStateToAmanzi(col,
                            bgc_props_, bgc_state_, bgc_aux_data_);

  return num_iterations;
}

} // namespace EcoSIM
} // namespace Amanzi
