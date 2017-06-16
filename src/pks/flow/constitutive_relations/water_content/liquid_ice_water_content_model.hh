/*
  The liquid + ice water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Water content for a two-phase, liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_MODEL_HH_
#define AMANZI_FLOW_LIQUID_ICE_WATER_CONTENT_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class LiquidIceWaterContentModel {

 public:
  explicit
  LiquidIceWaterContentModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double nl, double si, double ni, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double nl, double si, double ni, double cv) const;
  double DWaterContentDSaturationLiquid(double phi, double sl, double nl, double si, double ni, double cv) const;
  double DWaterContentDMolarDensityLiquid(double phi, double sl, double nl, double si, double ni, double cv) const;
  double DWaterContentDSaturationIce(double phi, double sl, double nl, double si, double ni, double cv) const;
  double DWaterContentDMolarDensityIce(double phi, double sl, double nl, double si, double ni, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double nl, double si, double ni, double cv) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif