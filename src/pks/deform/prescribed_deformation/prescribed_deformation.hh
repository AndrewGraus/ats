/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Prescribed Deformation PK.

   <ParameterList name="prescribed deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#ifndef PKS_PRESCRIBED_DEFORMATION_HH_
#define PKS_PRESCRIBED_DEFORMATION_HH_

#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "pk_physical_base.hh"

namespace Amanzi {
namespace Deform {

class PrescribedDeformation : public PKPhysicalBase {

 public:

  PrescribedDeformation(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                        Teuchos::ParameterList& FElist,
                        const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PrescribedDeformation() {}

  // ConstantTemperature is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods
  virtual bool advance(double dt);

  virtual double get_dt() {
    return 1.0e99;
  }

 private:

  int prescribed_deformation_case_;
  double deformation_fn_2(double x, double y, double t);
  double deformation_fn_3(double x, double y, double z, double t);
  double deformation_fn_4(double x, double y, double z, double t);
  double deformation_fn_5(double x, double y, double z, double t);
  double tmax_, sigma_, mag_, z0_, maxpert_;

  Key poro_key_;

  // A few options for advance
  bool advance_analytic_(double dt);

  // misc setup information
  Teuchos::ParameterList prescribed_deformation_plist_;

  // factory registration
  static RegisteredPKFactory<PrescribedDeformation> reg_;
};

} // namespace
} // namespace

#endif
