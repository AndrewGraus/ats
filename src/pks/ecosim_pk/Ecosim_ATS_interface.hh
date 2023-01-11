/* -*-  mode: c++; indent-tabs-mode: nil -*- */

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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include "BGCEngine.hh"
#include "PK_Factory.hh"
#include "pk_physical_default.hh"

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
  virtual ~EcoSIM() override {}

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

  void CopyToEcoSIM(int col,
          BGCProperties& props,
          BGCState& state,
          BGCAuxiliaryData& aux_data);

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

   int InitializeSingleColumn(int col, const std::string& condition);

   int AdvanceSingleColumn(double dt, int col);

   void CopyEcoSIMStateToAmanzi(
       const int cell,
       const BGCProperties& props,
       const BGCState& state,
       const BGCAuxiliaryData& aux_data);

   void ComputeNextTimeStep();

 protected:
  double dt_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_; //might need this?
  //Key domain_surf_;

  //The helper functions from BGC are protected not private (unclear why)
  //I don't think I need this here, probably in the engine
  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
          Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy);

  void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                                         double* col_vec, int ncol);
  void EcoSIM::ColDepthDz_(AmanziMesh::Entity_ID col,
                              Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                              Teuchos::Ptr<Epetra_SerialDenseVector> dz)

  //evaluator for transpiration;
  //I don't think I need this anymore
  //Teuchos::RCP<PrimaryVariableFieldEvaluator> p_root_eval_;

  // keys
  Key wc_key_;
  Key wc_root_key_;
  Key trans_key_;
  Key p_root_key_;
 private:
  //factory registration
  static RegisteredPKFactory<EcoSIM> reg_;
};

} // namespace EcoSIM
} // namespace Amanzi

#endif
