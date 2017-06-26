/*
  The three phase water content model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Water content for a three-phase, gas+liquid+ice evaluator.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_MODEL_HH_
#define AMANZI_FLOW_THREE_PHASE_WATER_CONTENT_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class ThreePhaseWaterContentModel {

 public:
  explicit
  ThreePhaseWaterContentModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDSaturationLiquid(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDMolarDensityLiquid(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDSaturationIce(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDMolarDensityIce(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDSaturationGas(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDMolarDensityGas(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDMolFracGas(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double nl, double si, double ni, double sg, double ng, double omega, double cv) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif