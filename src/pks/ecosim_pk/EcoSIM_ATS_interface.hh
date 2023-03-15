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

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "Key.hh"
#include "Mesh.hh"
#include "State.hh"
#include "BGCEngine.hh"
#include "PK_Factory.hh"
#include "pk_physical_default.hh"
#include "PK_Physical.hh"


namespace Amanzi {
namespace EcoSIM {

class EcoSIM : public PK_Physical_Default {

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
  virtual void Setup() override;

  // -- initalize owned (dependent) variables
  //virtual void Initialize(const Teuchos::Ptr<State>& S);
  virtual void Initialize() override;

  // --provide timestep size
  virtual double get_dt() override {
    return dt_;
  }

  virtual void set_dt(double dt) override {
    dt_ = dt;
  }

  // -- commit the model
  //virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- Update diagnostics for vis.
  //virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}
  //virtual void CalculateDiagnostics(const Tag& tag) override;

  // -- advance the model
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

  virtual std::string name(){return "EcoSIM for ATS";};

  /*void CopyToEcoSIM(int col,
          BGCProperties& props,
          BGCState& state,
          BGCAuxiliaryData& aux_data);*/

 private:

   //Helper functions from Alquimia
   void CopyToEcoSIM(int col,
           BGCProperties& props,
           BGCState& state,
           BGCAuxiliaryData& aux_data);

   void CopyFromEcoSIM(const int cell,
                const BGCProperties& props,
                const BGCState& state,
                const BGCAuxiliaryData& aux_data);

   int InitializeSingleColumn(int col);

   int AdvanceSingleColumn(double dt, int col);

   void CopyEcoSIMStateToAmanzi(
       const int cell,
       const BGCProperties& props,
       const BGCState& state,
       const BGCAuxiliaryData& aux_data);

   void ComputeNextTimeStep();

 protected:
  double dt_;
  double c_m_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_; //might need this?
  Key domain_surf_;

  //The helper functions from BGC are protected not private (unclear why)
  //I don't think I need this here, probably in the engine
  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy);

  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                                         double* col_vec, int ncol);
  void ColDepthDz_(AmanziMesh::Entity_ID col,
                              Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                              Teuchos::Ptr<Epetra_SerialDenseVector> dz);

  void ColumnToField_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                                 double* col_vec, int ncol);
  //evaluator for transpiration;
  //I don't think I need this anymore
  //Teuchos::RCP<PrimaryVariableFieldEvaluator> p_root_eval_;

  int number_aqueous_components_;
  int ncells_per_col_;
  int num_cols_;
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
  Key rock_den_key_;
  Key T_key_;
  Key conductivity_key_;
  Key cv_key_;
  Key min_vol_frac_key_;
  Key ecosim_aux_data_key_;

 private:
  BGCState bgc_state_;
  BGCProperties bgc_props_;
  BGCAuxiliaryData bgc_aux_data_;

 private:
  //factory registration
  static RegisteredPKFactory<EcoSIM> reg_;
};

} // namespace EcoSIM
} // namespace Amanzi

#endif
