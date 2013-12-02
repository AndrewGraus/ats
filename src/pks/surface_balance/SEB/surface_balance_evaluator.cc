/*
  SurfaceBalanceEvaluator evaluates SEB as a nonlinear source instead of as a PK.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Debugger.hh"
#include "surface_balance_evaluator.hh"
#include "SnowEnergyBalance.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceEvaluator::SurfaceBalanceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // my keys
  my_key_ = "surface_conducted_energy_source";

  // parameters
  min_wind_speed_ = plist_.get<double>("minimum wind speed", 1.0);
  snow_ground_trans_ = plist_.get<double>("minimum snow depth", 0.02);
  albedo_trans_ = plist_.get<double>("albedo transition depth", 0.02);

  // dependencies
  dependencies_.insert("surface_temperature");
  dependencies_.insert("snow_depth");
  dependencies_.insert("ponded_depth");
  dependencies_.insert("surface_porosity");
  dependencies_.insert("surface_vapor_pressure");
  dependencies_.insert("air_temperature");
  dependencies_.insert("incoming_shortwave_radiation");
  dependencies_.insert("relative_humidity");
  dependencies_.insert("wind_speed");
  dependencies_.insert("precipitation_rain");
  dependencies_.insert("precipitation_snow");

  //  dependencies_.insert("snow_density");
  //  dependencies_.insert("days_of_nosnow");
  //  dependencies_.insert("snow_temperature");

};


SurfaceBalanceEvaluator::SurfaceBalanceEvaluator(const SurfaceBalanceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    min_wind_speed_(other.min_wind_speed_),
    snow_ground_trans_(other.snow_ground_trans_),
    albedo_trans_(other.albedo_trans_),
    db_(other.db_) {}


Teuchos::RCP<FieldEvaluator>
SurfaceBalanceEvaluator::Clone() const {
  return Teuchos::rcp(new SurfaceBalanceEvaluator(*this));
}


void
SurfaceBalanceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  if (db_ == Teuchos::null) {
    // Debugger
    db_ = Teuchos::rcp(new Debugger(S->GetMesh("surface"), my_key_, plist_));
  }

  // Pull dependencies out of state.
  const Epetra_MultiVector& snow_temp = *S->GetFieldData("snow_temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_depth = *S->GetFieldData("snow_depth")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_density = *S->GetFieldData("snow_density")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& days_of_nosnow = *S->GetFieldData("days_of_nosnow")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_temp = *S->GetFieldData("surface_temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& ponded_depth = *S->GetFieldData("ponded_depth")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& surf_porosity = *S->GetFieldData("surface_porosity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& soil_vapor_pressure = *S->GetFieldData("surface_vapor_pressure")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& air_temp = *S->GetFieldData("air_temperature")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& incoming_shortwave = *S->GetFieldData("incoming_shortwave_radiation")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& relative_humidity = *S->GetFieldData("relative_humidity")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& wind_speed = *S->GetFieldData("wind_speed")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& precip_rain = *S->GetFieldData("precipitation_rain")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& precip_snow = *S->GetFieldData("precipitation_snow")
      ->ViewComponent("cell",false);

  Epetra_MultiVector& Qe = *result->ViewComponent("cell",false);
  
  // Create the SEB data structure
  SurfaceEnergyBalance::LocalData data;
  data.st_energy.dt = 0.;
  data.st_energy.AlbedoTrans = albedo_trans_;

  SurfaceEnergyBalance::LocalData data_bare;
  data_bare.st_energy.dt = 0.;
  data_bare.st_energy.AlbedoTrans = albedo_trans_;

  
  int count = Qe.MyLength();
  for (unsigned int c=0; c!=count; ++c) {
    // ATS Calcualted Data
    double density_air = 1.275;       // Density of Air ------------------- [kg/m^3]
    data.st_energy.water_depth = ponded_depth[0][c];
    data.st_energy.temp_ground = surf_temp[0][c];
    data.vp_ground.temp = surf_temp[0][c];
    data.vp_ground.actual_vaporpressure = soil_vapor_pressure[0][c] * data.st_energy.Apa; //Converts Mole fraction to vapor pressure [moleFraction/atmosphericPressure] 
    data.st_energy.porrowaLe = surf_porosity[0][c] * density_air * data.st_energy.Le;
    // MET station data
    data.st_energy.temp_air = air_temp[0][c];
    data.st_energy.QswIn = incoming_shortwave[0][c];
    data.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
    data.vp_air.temp = air_temp[0][c];
    data.vp_air.relative_humidity = relative_humidity[0][c];
    // STORED INFO FOR SnowEnergyBalanc Model
    data.st_energy.ht_snow = snow_depth[0][c];
    data.st_energy.density_snow = snow_density[0][c];
    data.st_energy.age_snow = days_of_nosnow[0][c];

    // Extras just for the evaluator
    data.st_energy.temp_snow = snow_temp[0][c];

    std::cout << "Updating SEB: ht_snow, tmp_snow = " << data.st_energy.ht_snow << ", " << data.st_energy.temp_snow << std::endl;

    // Snow-ground Smoothing 
    if ((data.st_energy.ht_snow > snow_ground_trans_) ||
        (data.st_energy.ht_snow <= 0)) {
      // Run the Snow Energy Balance Model as normal.
      SurfaceEnergyBalance::UpdateEnergyBalance(data);

      std::cout << "  SEB, non-averaged: ht_snow, tmp_snow = " << data.st_energy.ht_snow << ", " << data.st_energy.temp_snow << std::endl;

    } else {
      double theta=0.0;

      // Fraction of snow covered ground
      //      theta = pow ((data.st_energy.ht_snow / snow_ground_trans_),2);
      theta = data.st_energy.ht_snow / snow_ground_trans_;

      // Calculate as if ht_snow is the min value.
      data.st_energy.ht_snow = snow_ground_trans_;
      SurfaceEnergyBalance::UpdateEnergyBalance(data);

      // Calculate as if bare ground
      // ATS Calcualted Data
      data_bare.st_energy.water_depth = ponded_depth[0][c];
      data_bare.st_energy.temp_ground = surf_temp[0][c];
      data_bare.vp_ground.temp = surf_temp[0][c];
      data_bare.vp_ground.actual_vaporpressure = soil_vapor_pressure[0][c] * data_bare.st_energy.Apa; //Converts Mole fraction to vapor pressure [moleFraction/atmosphericPressure]
      data_bare.st_energy.porrowaLe = surf_porosity[0][c] * density_air * data_bare.st_energy.Le;
      // MET station data
      data_bare.st_energy.temp_air = air_temp[0][c];
      data_bare.st_energy.QswIn = incoming_shortwave[0][c];
      data_bare.st_energy.Us = std::max(wind_speed[0][c], min_wind_speed_);
      data_bare.vp_air.temp = air_temp[0][c];
      data_bare.vp_air.relative_humidity = relative_humidity[0][c];
      // STORED INFO FOR SnowEnergyBalanc Model
      data_bare.st_energy.ht_snow = snow_depth[0][c];
      data_bare.st_energy.density_snow = snow_density[0][c];
      data_bare.st_energy.age_snow = days_of_nosnow[0][c];
      data_bare.st_energy.ht_snow = 0.;
      SurfaceEnergyBalance::UpdateEnergyBalance(data_bare);

      // Calculating Data for ATS
      std::cout << "  SEB, averaged: ht_snow, tmp_snow = " << data.st_energy.ht_snow << ", " << data.st_energy.temp_snow << std::endl;
      std::cout << "  Averaging fQc (theta=" << theta << "): " << data.st_energy.fQc << ", " << data_bare.st_energy.fQc << std::endl;
      data.st_energy.fQc = data.st_energy.fQc * theta + data_bare.st_energy.fQc * (1.-theta);
    }

    // Store results
    Qe[0][c] = data.st_energy.fQc;
  }

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("air_temp"); vecs.push_back(S->GetFieldData("air_temperature").ptr());
    vnames.push_back("air_temp2"); vecs.push_back(S->GetFieldData("air_temperature").ptr());
    vnames.push_back("Qsw_in"); vecs.push_back(S->GetFieldData("incoming_shortwave_radiation").ptr());
    vnames.push_back("precip_rain"); vecs.push_back(S->GetFieldData("precipitation_rain").ptr());
    vnames.push_back("precip_snow"); vecs.push_back(S->GetFieldData("precipitation_snow").ptr());
    vnames.push_back("soil vapor pressure"); vecs.push_back(S->GetFieldData("surface_vapor_pressure").ptr());
    vnames.push_back("T_ground"); vecs.push_back(S->GetFieldData("surface_temperature").ptr());
    vnames.push_back("water_source"); vecs.push_back(S->GetFieldData("surface_mass_source").ptr());
    vnames.push_back("e_source"); vecs.push_back(S->GetFieldData("surface_conducted_energy_source").ptr());
    db_->WriteVectors(vnames, vecs, false);
  }
}


void SurfaceBalanceEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& results) {
  results->PutScalar(0.);
}

} // namespace
} // namespac
