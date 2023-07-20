/*--------------------------------------------------------------------------
  ATS

  License: see $ATS_DIR/COPYRIGHT
  Author: Andrew Graus

  The idea here is to begin writing the EcoSIM ATS interface with a simple
  program. To start we are going to try to do a few things:

  1) Initalize a PK called EcoSIM_ATS
  2) Have that PK take in the water content
  3) modify the water content in a simple way to mock roots (take away water)
  4) modify it so it will take in tracers (how roots take in nutrients)

  --------------------------------------------------------------------------*/
//Eventually add if statement here (probably tied to something at compile time
//
//#ifndef PKS_BGC_SIMPLE_HH_
//#define PKS_BGC_SIMPLE_HH_

#ifndef PKS_ECOSIM_HH_
#define PKS_ECOSIM_HH_

#include <map>
#include <vector>
#include <string>

#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "Key.hh"
#include "Mesh.hh"
#include "State.hh"
#include "BGCEngine.hh"
#include "PK_Factory.hh"
#include "pk_physical_default.hh"
#include "PK_Physical.hh"
#include "ecosim_mod_test_wrapper.h"


namespace Amanzi {
namespace EcoSIM {

class EcoSIM : public PK_Physical {

 public:

  //Unclear if the constructor is neccessary
  EcoSIM(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& plist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  ~EcoSIM();

  // is a PK
  // -- Setup data
  //virtual void Setup(const Teuchos::Ptr<State>&S);
  virtual void Setup() final;

  // -- initalize owned (dependent) variables
  //virtual void Initialize(const Teuchos::Ptr<State>& S);
  virtual void Initialize() final;

  // --provide timestep size
  virtual double get_dt() final {
    return dt_;
  }

  virtual void set_dt(double dt) final {
    dt_ = dt;
  }

  // -- commit the model
  //virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;

  // -- Update diagnostics for vis.
  //virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}
  //virtual void CalculateDiagnostics(const Tag& tag) override;

  // -- advance the model
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) final;

  virtual std::string name(){return "EcoSIM for ATS";};

  //This is not in the Alquimia_PK, for whatever reason it is defined in
  //The Chemistry_PK even though it isn't used there, and then included
  Teuchos::RCP<BGCEngine> bgc_engine() { return bgc_engine_; }

  /*void CopyToEcoSIM(int col,
          BGCProperties& props,
          BGCState& state,
          BGCAuxiliaryData& aux_data);*/

 private:

   //Helper functions from Alquimia
   void CopyToEcoSIM(int col,
           BGCProperties& props,
           BGCState& state,
           BGCAuxiliaryData& aux_data,
         const Tag& water_tag = Tags::DEFAULT);

   void CopyFromEcoSIM(const int cell,
                const BGCProperties& props,
                const BGCState& state,
                const BGCAuxiliaryData& aux_data,
              const Tag& water_tag = Tags::DEFAULT);

    //Helper functions from Alquimia
    void CopyToEcoSIM_process(int proc,
            BGCProperties& props,
            BGCState& state,
            BGCAuxiliaryData& aux_data,
          const Tag& water_tag = Tags::DEFAULT);

    void CopyFromEcoSIM_process(const int proc,
                 const BGCProperties& props,
                 const BGCState& state,
                 const BGCAuxiliaryData& aux_data,
               const Tag& water_tag = Tags::DEFAULT);

   int InitializeSingleColumn(int col);

   int AdvanceSingleColumn(double dt, int col);

   int InitializeSingleProcess(int proc);

   int AdvanceSingleProcess(double dt, int proc);

   void CopyEcoSIMStateToAmanzi(
       const int cell,
       const BGCProperties& props,
       const BGCState& state,
       const BGCAuxiliaryData& aux_data,
     const Tag& water_tag = Tags::DEFAULT);

   void ComputeNextTimeStep();

 protected:
  double dt_;
  double c_m_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_; //might need this?
  Key domain_surf_;
  std::string passwd_ = "state";

  //The helper functions from BGC are protected not private (unclear why)
  //I don't think I need this here, probably in the engine
  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void FieldToColumn_(AmanziMesh::Entity_ID col, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void MatrixFieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_MultiVector& m_arr,
      Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr);

  //void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
  //                                       double* col_vec);
  void ColDepthDz_(AmanziMesh::Entity_ID col,
                              Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                              Teuchos::Ptr<Epetra_SerialDenseVector> dz);

  void ColumnToField_(AmanziMesh::Entity_ID col, Epetra_Vector& vec,
                                 Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void ColumnToField_(AmanziMesh::Entity_ID col, Teuchos::Ptr<Epetra_SerialDenseVector> vec,
                                 Teuchos::Ptr<Epetra_SerialDenseVector> col_vec);

  void MatrixColumnToField_(AmanziMesh::Entity_ID col, Epetra_MultiVector& m_arr,
                                 Teuchos::Ptr<Epetra_SerialDenseMatrix> col_arr);
  //evaluator for transpiration;
  //I don't think I need this anymore
  //Teuchos::RCP<PrimaryVariableFieldEvaluator> p_root_eval_;

  int number_aqueous_components_;
  int ncells_per_col_;
  int num_cols_;
  int ncols_local;
  int ncols_global;
  int ncols_global_ptype;
  double saved_time_;
  double current_time_;

  // keys
  Key tcc_key_;
  Key poro_key_;
  Key saturation_liquid_key_;
  Key saturation_gas_key_;
  Key saturation_ice_key_;
  Key elev_key_;
  Key water_content_key_;
  Key rel_perm_key_;
  Key liquid_den_key_;
  Key ice_den_key_;
  Key gas_den_key_;
  Key gas_den_key_test_;
  Key rock_den_key_;
  Key T_key_;
  Key therm_cond_key_;
  Key cv_key_;
  Key min_vol_frac_key_;
  Key ecosim_aux_data_key_;
  Key bulk_dens_key_;
  Key hydra_cond_key_;
  Key sw_key_;
  Key lw_key_;
  Key air_temp_key_;
  Key vp_air_key_;
  Key wind_speed_key_;
  Key prain_key_;
  Key f_wp_key_;
  Key f_root_key_;
  Key suc_key_;
  Key aspect_key_;
  Key slope_key_;

  Teuchos::RCP<BGCEngine> bgc_engine_;

  double atm_n2_, atm_o2_, atm_co2_, atm_ch4_, atm_n2o_, atm_h2_, atm_nh3_;

 private:
  BGCState bgc_state_;
  BGCProperties bgc_props_;
  BGCAuxiliaryData bgc_aux_data_;
  BGCSizes bgc_sizes_;

  Teuchos::RCP<Epetra_SerialDenseVector> col_vol_save;
  Teuchos::RCP<Epetra_SerialDenseVector> col_wc_save;

  bool bgc_initialized_;
  bool has_energy, has_gas, has_ice;
  std::vector<std::string> component_names_;
  int num_components;

 private:
  //factory registration
  static RegisteredPKFactory<EcoSIM> reg_;
};

} // namespace EcoSIM
} // namespace Amanzi

#endif
