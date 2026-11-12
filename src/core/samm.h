/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
SAMM - Soil Aggregation through Microbial Mediation, version 1.1
Laub, M. et al. (2024), Geosci. Model Dev. 17, 931-956, https://doi.org/10.5194/gmd-17-931-2024,
corrigendum (2025) and R code v1.1 (https://doi.org/10.5281/zenodo.15648425).

Mineral N coupling (gross mineralisation of microbial N surplus and turnover to mineral N, optional
microbial immobilisation of mineral N) and MAOM desorption as in the "Samm3H_functionsMinN.R"
scenario D of the SAMM ecosystem-model development.

This file has no MONICA dependencies, so the same core can be tested against the R model.
All pools are in kg ha-1 per layer, rates per day.
*/

#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace samm {

enum Pool : int {
  STR_C = 0,  //!< structural litter C (lignin + polyphenols)
  LAB_C,      //!< labile (metabolic) litter C
  LMW_C,      //!< low molecular weight (dissolved) C
  MIC_C,      //!< microbial biomass C
  MAO_C,      //!< mineral-associated C
  AggSTR_C,   //!< structural litter C protected in aggregates
  AggLAB_C,   //!< labile litter C protected in aggregates
  AggMAO_C,   //!< mineral-associated C protected in aggregates
  LAB_N,
  LMW_N,
  MIC_N,
  MAO_N,
  AggLAB_N,
  AggMAO_N,
  NUM_POOLS
};

const char* poolName(int pool);

using LayerPools = std::array<double, NUM_POOLS>;

//! cumulative (or daily) fluxes of one layer, kg ha-1
enum Flux : int {
  F_CO2 = 0,        //!< heterotrophic CO2-C
  F_LEACH_C_OUT,    //!< LMW-C moved across the lower boundary of the layer (positive downward)
  F_LEACH_N_OUT,    //!< LMW-N moved across the lower boundary of the layer (positive downward)
  F_MIN_GROSS,      //!< N moved from LMW-N to mineral N (gross mineralisation)
  F_MIN_IMM,        //!< mineral N immobilised by microbes
  F_INPUT_C,        //!< C added through input rates (only used by rate-based inputs)
  F_INPUT_N,        //!< N added through input rates
  NUM_FLUXES
};

using LayerFluxes = std::array<double, NUM_FLUXES>;

struct Params {
  // turnover rates [d-1] at rate modifier 1
  double k_STR{0.0031};
  double k_LAB{0.0366};
  double k_MIC{0.007};
  double k_MAO{0.00041};
  double mu_max{0.337};
  double kAgg{0.0323};
  double K_M_MIC_C{66.8};     //!< half saturation of the microbial activity factor [kg C ha-1 per refDepth]
  double m_MIC{0.00064};
  double K_LMW_MAO{0.018};
  double c_SORP{0.049};       //!< maximum sorption capacity [kg C kg-1 silt+clay]
  // carbon use efficiencies
  double CUE_STR{0.57};
  double CUE_LAB{0.51};
  double CUE_LMW{0.6};
  // other parameters
  double MicMinCN{5.72};
  double MicMaxCN{7.77};
  double f_MICMAOM{0.21};
  double pc_STR_LAB{2.54};
  double aggfactSTR_C{1.16};
  double aggfactMAO_C{1.57};
  double NonMicAgg{32.5};      //!< physicochemical aggregate formation [kg MIC-C eq ha-1 d-1 per refDepth]
  double refDepth{0.15};       //!< depth K_M_MIC_C and NonMicAgg refer to [m]
  //! if > 0: physicochemical aggregate formation follows the rate modifier, NonMicAgg * corr_F / NonMicAggCorrFRef,
  //! like aggregate breakdown; the value is the mean corr_F at which NonMicAgg was calibrated (1.942 for Khon Kaen).
  //! 0 (default) = constant formation as in SAMM v1.1
  double NonMicAggCorrFRef{0.0};
  double a_MIC_min{0.05};
  double leachCap{0.95};       //!< maximum share of a donor pool moved by water per day
  //! adaptive sub-steps (Profile::stepAdaptive): the day is split into n steps so that no pool loses more than this
  //! share per step at the start of the day; a step whose RK4 stages would make a pool negative is redone with twice
  //! as many steps. 0 = always one RK4 step per day (exactly as deSolve::rk4 in the R code)
  double maxTurnoverPerStep{1.0};
  double maxSubsteps{64};      //!< upper limit of the number of steps per day
  // mineral N coupling and desorption (scenario D)
  double f_min{0.5};           //!< share of microbial N surplus going to mineral N
  double f_turn{0.10};         //!< share of MIC-N and MAO-N turnover to LMW-N going to mineral N
  double min_imm_frac{0.0};    //!< daily share of mineral N microbes may take when N-limited
  double min_imm_k{0.0};       //!< first-order microbial uptake of mineral N [d-1]
  double k_des{0.006};         //!< MAOM desorption [kg C t-1 soil d-1] at full saturation

  //! SAMM v1.1 with aggregates (corrigendum Table 4, SAMM column), scenario D mineral N coupling
  static Params samm_v11();
  //! SAMMnoAgg v1.1 (corrigendum Table 4, SAMMnoAgg column), scenario D mineral N coupling
  static Params sammNoAgg_v11();
  //! switch off the scenario D extensions (pure SAMM v1.1 as in the R code)
  Params& withoutMineralN() { f_min = 0; f_turn = 0; min_imm_frac = 0; min_imm_k = 0; k_des = 0; return *this; }
};

//! static properties of a layer
struct LayerEnv {
  double thickness{0.1};      //!< [m]
  double bulkDensity{1450};   //!< [kg m-3]
  double siltClay{0.2};       //!< [kg kg-1]
  double cSorp{-1};           //!< layer specific c_SORP [kg C kg-1 silt+clay], < 0: use Params::c_SORP

  double soilMass_t_ha() const { return thickness * bulkDensity / 1000.0 * 10000.0; }
  //! maximum MAO-C (corrected equation of the corrigendum) [kg C ha-1]
  double maoCmax(const Params& p) const {
    return thickness * 10000.0 * bulkDensity * siltClay * (cSorp >= 0 ? cSorp : p.c_SORP);
  }
};

//! drivers of one layer for one evaluation of the derivatives
struct LayerDrivers {
  double corrF{1.0};          //!< combined temperature x moisture rate modifier
  double water_kg_ha{1.0e5};  //!< soil water in the layer [kg ha-1]
  double minN_kg_ha{0.0};     //!< mineral N available for microbial immobilisation [kg N ha-1]
};

struct Drivers {
  std::vector<LayerDrivers> layers;
  //! water flux at the lower boundary of each layer [mm d-1], positive downward
  std::vector<double> waterFlux_mm;
  double deepConcC{0.0};      //!< LMW-C concentration of water entering from below [kg kg-1]
  double deepConcN{0.0};
  //! optional continuous inputs [kg ha-1 d-1] per layer and pool (empty = none)
  std::vector<LayerPools> inputRates;
};

class Profile {
public:
  Profile() = default;
  Profile(Params p, std::vector<LayerEnv> env) : _p(p), _env(std::move(env)) {}

  const Params& params() const { return _p; }
  Params& params() { return _p; }
  const std::vector<LayerEnv>& env() const { return _env; }
  std::vector<LayerEnv>& env() { return _env; }
  size_t size() const { return _env.size(); }

  //! derivatives of all pools and flux rates
  void derivatives(const std::vector<LayerPools>& y, const Drivers& d,
                   std::vector<LayerPools>& dy, std::vector<LayerFluxes>& fluxRates) const;

  //! one classical 4th order Runge-Kutta step of dt days starting at t0 (days since the start of the day);
  //! dt = 1 is deSolve::rk4 with daily times. drivers(t) is called with t0, t0 + dt/2 (twice) and t0 + dt.
  //! The fluxes of the step are added to fluxes. If negativeStage is given, it is set to true when a stage state
  //! or the result has a negative pool.
  void stepRK4(std::vector<LayerPools>& y, const std::function<const Drivers&(double)>& drivers,
               std::vector<LayerFluxes>& fluxes, double t0 = 0.0, double dt = 1.0, bool* negativeStage = nullptr) const;

  //! same, with drivers constant over the day
  void stepRK4(std::vector<LayerPools>& y, const Drivers& d, std::vector<LayerFluxes>& fluxes) const;

  //! one day in n equal RK4 steps, n chosen as described for Params::maxTurnoverPerStep; returns n
  int stepAdaptive(std::vector<LayerPools>& y, const std::function<const Drivers&(double)>& drivers,
                   std::vector<LayerFluxes>& fluxes) const;
  int stepAdaptive(std::vector<LayerPools>& y, const Drivers& d, std::vector<LayerFluxes>& fluxes) const;

  //! carbon and nitrogen balance residual of one step [kg ha-1] (should be ~0)
  static std::pair<double, double> balanceResidual(const std::vector<LayerPools>& before,
                                                   const std::vector<LayerPools>& after,
                                                   const std::vector<LayerFluxes>& dailyFluxes);
private:
  Params _p;
  std::vector<LayerEnv> _env;
};

double totalC(const LayerPools& y);
double totalN(const LayerPools& y);
//! C that MONICA reports as soil organic carbon: all pools except free (particulate) litter STR_C and LAB_C
double socC(const LayerPools& y);

//! rate modifiers of the SAMM paper (reference 20 C)
double rateModTempSAMM(double temperature);
double rateModWaterSAMM(double pF);
//! temperature function of Kirschbaum (1995, Soil Biol. Biochem. 27, 753-760), maximum 1 at 36.9 C
double rateModTempKirschbaum(double temperature);

} // namespace samm
