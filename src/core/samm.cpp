/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
SAMM v1.1 core, see samm.h. The equations follow Samm_main.R (SAMM v1.1) and, for the mineral N
coupling and desorption, Samm3H_functionsMinN.R. Variable names follow the R code.
*/

#include "samm.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace samm {

const char* poolName(int pool) {
  static const char* names[NUM_POOLS] = {"STR_C", "LAB_C", "LMW_C", "MIC_C", "MAO_C", "AggSTR_C", "AggLAB_C",
                                          "AggMAO_C", "LAB_N", "LMW_N", "MIC_N", "MAO_N", "AggLAB_N", "AggMAO_N"};
  return pool >= 0 && pool < NUM_POOLS ? names[pool] : "";
}

Params Params::samm_v11() {
  return Params();
}

Params Params::sammNoAgg_v11() {
  Params p;
  p.k_STR = 0.0032;
  p.k_LAB = 0.0378;
  p.k_MIC = 0.0068;
  p.k_MAO = 0.00028;
  p.mu_max = 0.29;
  p.kAgg = 1.0;           // set to 1 to deactivate aggregates
  p.K_M_MIC_C = 49.6;
  p.m_MIC = 0.0023;
  p.K_LMW_MAO = 0.067;
  p.CUE_STR = 0.66;
  p.CUE_LAB = 0.72;
  p.MicMinCN = 5.63;
  p.MicMaxCN = 7.31;
  p.f_MICMAOM = 0.17;
  p.pc_STR_LAB = 3.94;
  p.aggfactSTR_C = 0.0;   // set to 0 to deactivate aggregates
  p.aggfactMAO_C = 0.0;
  p.NonMicAgg = 0.0;
  return p;
}

double totalC(const LayerPools& y) {
  return y[STR_C] + y[LAB_C] + y[LMW_C] + y[MIC_C] + y[MAO_C] + y[AggSTR_C] + y[AggLAB_C] + y[AggMAO_C];
}

double totalN(const LayerPools& y) {
  return y[LAB_N] + y[LMW_N] + y[MIC_N] + y[MAO_N] + y[AggLAB_N] + y[AggMAO_N];
}

double socC(const LayerPools& y) {
  return y[LMW_C] + y[MIC_C] + y[MAO_C] + y[AggSTR_C] + y[AggLAB_C] + y[AggMAO_C];
}

double rateModTempSAMM(double temperature) {
  return std::pow(2.0, (temperature - 20.0) / 10.0);
}

double rateModTempKirschbaum(double temperature) {
  return std::exp(-3.764 + 0.204 * temperature * (1.0 - 0.5 * temperature / 36.9));
}

double rateModWaterSAMM(double pF) {
  double corr_M = 1;
  if (pF >= 2.5) corr_M = 1.625 - pF / 4.0;
  if (pF >= 6.5) corr_M = 0;
  if (pF <= 1.5) corr_M = 0.6 + 0.4 * pF / 1.5;
  if (pF <= 0) corr_M = 0.6;
  return corr_M;
}

namespace {

//! dynamic, C:N dependent carbon use efficiency of LMW (CUE_CN_LMW in the R code)
double cueCnLmw(double LMW_C, double LMW_N, double CUE_LMW) {
  // safety net: no growth from an empty (or, in an RK4 stage, negative) LMW pool (the R code has no lower limit)
  if (LMW_C <= 0 || LMW_N <= 0) return 0.0;
  double cue = 1.0 / ((1e-6 + LMW_C) / (LMW_N + 1e-10)) * 13.4;
  if (std::isnan(cue)) cue = 0;
  if (cue > 1) cue = 1;
  return CUE_LMW * cue;
}

double ratio(double a, double b) { return b != 0 ? a / b : 0.0; }

//! solute moved across a layer boundary by a signed water flux (LayerSoluteTransport in the R code),
//! positive downward
double soluteTransport(double waterFlux_mm, double upperPool, double upperWater,
                       const double* lowerPool, double lowerWater, double deepConc, double cap) {
  double flux_kg_ha = waterFlux_mm * 10000.0;
  if (flux_kg_ha >= 0) {
    double donor = std::max(upperPool, 0.0);
    double water = std::max(upperWater, std::numeric_limits<double>::epsilon());
    return std::min(flux_kg_ha * donor / water, cap * donor);
  }
  if (!lowerPool) return flux_kg_ha * std::max(deepConc, 0.0);
  double donor = std::max(*lowerPool, 0.0);
  double water = std::max(lowerWater, std::numeric_limits<double>::epsilon());
  return -std::min(std::abs(flux_kg_ha) * donor / water, cap * donor);
}

} // namespace

void Profile::derivatives(const std::vector<LayerPools>& y, const Drivers& drv,
                          std::vector<LayerPools>& dy, std::vector<LayerFluxes>& fr) const {
  const auto& p = _p;
  const size_t n = _env.size();
  dy.assign(n, LayerPools{});
  fr.assign(n, LayerFluxes{});

  for (size_t i = 0; i < n; i++) {
    const auto& s = y[i];
    auto& d = dy[i];
    const auto& env = _env[i];
    const auto& ld = drv.layers.at(i);
    const double corr_F = ld.corrF;
    const double depthScale = env.thickness / p.refDepth;
    const double K_M_MIC_C = p.K_M_MIC_C * depthScale;
    const double NonMicAgg = p.NonMicAgg * depthScale * (p.NonMicAggCorrFRef > 0 ? corr_F / p.NonMicAggCorrFRef : 1.0);
    const double MAO_C_max = env.maoCmax(p);

    const double STR_Cv = s[STR_C], LAB_Cv = s[LAB_C], LMW_Cv = s[LMW_C], MIC_Cv = s[MIC_C], MAO_Cv = s[MAO_C];
    const double AggSTR_Cv = s[AggSTR_C], AggLAB_Cv = s[AggLAB_C], AggMAO_Cv = s[AggMAO_C];
    const double LAB_Nv = s[LAB_N], LMW_Nv = s[LMW_N], MIC_Nv = s[MIC_N], MAO_Nv = s[MAO_N];
    const double AggLAB_Nv = s[AggLAB_N], AggMAO_Nv = s[AggMAO_N];

    // helper variables
    double p_LAB = STR_Cv * p.pc_STR_LAB / std::max(LAB_Cv, 0.0000001);
    if (p_LAB > 1) p_LAB = 1;
    const double a_MIC = std::max(MIC_Cv / (K_M_MIC_C + MIC_Cv), p.a_MIC_min);
    const double cueLMW = cueCnLmw(LMW_Cv, LMW_Nv, p.CUE_LMW);

    // carbon flows
    const double STR_C_to_CO2 = STR_Cv * p.k_STR * a_MIC * (1 - p.CUE_STR) * corr_F;
    const double STR_C_to_LMW_C = STR_Cv * p.k_STR * a_MIC * (p.CUE_STR) * corr_F;
    const double LAB_C_to_CO2 = LAB_Cv * p.k_LAB * (1 - p_LAB) * a_MIC * (1 - p.CUE_LAB) * corr_F;
    const double LAB_C_to_LMW_C = LAB_Cv * p.k_LAB * a_MIC * (1 - p_LAB) * (p.CUE_LAB) * corr_F;
    const double LMW_C_to_CO2 = p.mu_max * LMW_Cv * a_MIC * (1 - cueLMW) * corr_F;
    const double LMW_C_to_MIC_C = p.mu_max * LMW_Cv * a_MIC * cueLMW * corr_F;
    // safety net: aggregate formation cannot become negative (only possible from negative RK4 stage states)
    const double STR_C_to_AggSTR_C = std::max(std::min((LMW_C_to_MIC_C + NonMicAgg) * p.aggfactSTR_C, STR_Cv), 0.0);
    const double LAB_C_to_AggLAB_C = std::max(std::min(STR_C_to_AggSTR_C * p.pc_STR_LAB, LAB_Cv), 0.0);
    const double MAO_C_to_AggMAO_C = std::max(std::min((LMW_C_to_MIC_C + NonMicAgg) * p.aggfactMAO_C, MAO_Cv), 0.0);
    const double AggSTR_C_to_STR_C = AggSTR_Cv * p.kAgg * corr_F;
    const double AggLAB_C_to_LAB_C = AggLAB_Cv * p.kAgg * corr_F;
    const double AggMAO_C_to_MAO_C = AggMAO_Cv * p.kAgg * corr_F;
    const double MIC_C_to_CO2 = MIC_Cv * p.m_MIC * corr_F;
    // microbial death, density dependent if betaMIC != 1 (rate k_MIC at MIC_C = MIC_ref)
    const double kDeath = p.betaMIC == 1.0 ? p.k_MIC
                          : p.k_MIC * std::pow(std::max(MIC_Cv, 0.0) / (p.MIC_ref * depthScale), p.betaMIC - 1.0);
    const double MIC_C_to_MAO_C = p.f_MICMAOM * MIC_Cv * kDeath * corr_F;
    const double MIC_C_to_LMW_C = (1 - p.f_MICMAOM) * MIC_Cv * kDeath * corr_F;
    const double MAO_C_to_LMW_C = MAO_Cv * p.k_MAO * a_MIC * corr_F;
    const double LMW_C_to_MAO_C = LMW_Cv * p.K_LMW_MAO * (MAO_C_max - MAO_Cv) / MAO_C_max * corr_F;
    // desorption as in COMISSION (Ahrens et al. 2015), scenario D
    const double MAO_C_to_LMW_CDesorption =
        p.k_des * env.soilMass_t_ha() * std::max(MAO_Cv, 0.0) / MAO_C_max * corr_F;

    d[STR_C] = -STR_C_to_LMW_C - STR_C_to_CO2 - STR_C_to_AggSTR_C + AggSTR_C_to_STR_C;
    d[LAB_C] = -LAB_C_to_LMW_C - LAB_C_to_CO2 - LAB_C_to_AggLAB_C + AggLAB_C_to_LAB_C;
    d[LMW_C] = +STR_C_to_LMW_C + LAB_C_to_LMW_C + MIC_C_to_LMW_C + MAO_C_to_LMW_C
               - LMW_C_to_CO2 - LMW_C_to_MIC_C - LMW_C_to_MAO_C + MAO_C_to_LMW_CDesorption;
    d[MIC_C] = +LMW_C_to_MIC_C - MIC_C_to_CO2 - MIC_C_to_MAO_C - MIC_C_to_LMW_C;
    d[MAO_C] = +MIC_C_to_MAO_C + LMW_C_to_MAO_C - MAO_C_to_LMW_CDesorption - MAO_C_to_LMW_C
               - MAO_C_to_AggMAO_C + AggMAO_C_to_MAO_C;
    d[AggSTR_C] = +STR_C_to_AggSTR_C - AggSTR_C_to_STR_C;
    d[AggLAB_C] = +LAB_C_to_AggLAB_C - AggLAB_C_to_LAB_C;
    d[AggMAO_C] = +MAO_C_to_AggMAO_C - AggMAO_C_to_MAO_C;
    fr[i][F_CO2] = MIC_C_to_CO2 + LMW_C_to_CO2 + LAB_C_to_CO2 + STR_C_to_CO2;

    // nitrogen flows
    double LAB_N_to_LMW_N = LAB_Nv * p.k_LAB * a_MIC * (1 - p_LAB) * corr_F;
    if (LAB_N_to_LMW_N < 0) LAB_N_to_LMW_N = 0;
    const double LMW_N_to_Mic_N = p.mu_max * LMW_Nv * a_MIC * corr_F;
    const double MIC_N_to_MAO_N = p.f_MICMAOM * MIC_Nv * kDeath * corr_F;
    const double MIC_N_to_LMW_N = (1 - p.f_MICMAOM) * MIC_Nv * kDeath * corr_F;
    const double MAO_N_to_LMW_N = MAO_Nv * p.k_MAO * a_MIC * corr_F;
    const double LMW_N_to_MAO_N = LMW_C_to_MAO_C * ratio(LMW_Nv, LMW_Cv);
    const double MAO_to_LMW_NDesorption = MAO_C_to_LMW_CDesorption * ratio(MAO_Nv, MAO_Cv);
    const double LAB_N_to_AggLAB_N = LAB_C_to_AggLAB_C * ratio(LAB_Nv, LAB_Cv);
    const double MAO_N_to_AggMAO_N = MAO_C_to_AggMAO_C * ratio(MAO_Nv, MAO_Cv);
    const double AggLAB_N_to_LAB_N = AggLAB_Nv * p.kAgg * corr_F;
    const double AggMAO_N_to_MAO_N = AggMAO_Nv * p.kAgg * corr_F;

    // LMW-N release and immobilization by microbes
    double MicLMW_NRelease = 0;
    const double MIC_CN = (MIC_Cv + 1E-10) / (MIC_Nv + 1E-11);
    if (MIC_CN < p.MicMinCN) MicLMW_NRelease = 0.5 * (MIC_Nv - (MIC_Cv / p.MicMinCN));
    double MicLMW_NImmobilization = 0;
    const double max_MicLMW_NImmobilization = LMW_Nv * 0.5;
    const double LMW_N_need = (MIC_Cv / p.MicMaxCN) - MIC_Nv;
    if (LMW_N_need > 0) MicLMW_NImmobilization = max_MicLMW_NImmobilization;
    if (max_MicLMW_NImmobilization > LMW_N_need) MicLMW_NImmobilization = LMW_N_need;

    // mineral N coupling (scenario D): gross mineralisation and immobilisation of mineral N
    const double minN = std::max(ld.minN_kg_ha, 0.0);
    const double Min_N_from_Mic = p.f_min * (MicLMW_NRelease + std::max(-MicLMW_NImmobilization, 0.0))
                                  + p.f_turn * (MIC_N_to_LMW_N + MAO_N_to_LMW_N);
    double Min_N_imm = 0;
    if (LMW_N_need > 0) {
      Min_N_imm = std::max(std::min(LMW_N_need - std::max(MicLMW_NImmobilization, 0.0), p.min_imm_frac * minN), 0.0);
    }
    Min_N_imm += std::max(std::min({p.min_imm_k * minN * a_MIC * corr_F, 0.5 * minN, MIC_Cv / p.MicMinCN - MIC_Nv}), 0.0);

    d[LAB_N] = -LAB_N_to_LMW_N - LAB_N_to_AggLAB_N + AggLAB_N_to_LAB_N;
    d[LMW_N] = +LAB_N_to_LMW_N - LMW_N_to_Mic_N + MIC_N_to_LMW_N + MAO_N_to_LMW_N + MicLMW_NRelease
               - MicLMW_NImmobilization - LMW_N_to_MAO_N + MAO_to_LMW_NDesorption - Min_N_from_Mic;
    d[MIC_N] = +LMW_N_to_Mic_N - MIC_N_to_MAO_N - MIC_N_to_LMW_N - MicLMW_NRelease + MicLMW_NImmobilization
               + Min_N_imm;
    d[MAO_N] = +MIC_N_to_MAO_N + LMW_N_to_MAO_N - MAO_to_LMW_NDesorption - MAO_N_to_LMW_N - MAO_N_to_AggMAO_N
               + AggMAO_N_to_MAO_N;
    d[AggLAB_N] = LAB_N_to_AggLAB_N - AggLAB_N_to_LAB_N;
    d[AggMAO_N] = MAO_N_to_AggMAO_N - AggMAO_N_to_MAO_N;
    fr[i][F_MIN_GROSS] = Min_N_from_Mic;
    fr[i][F_MIN_IMM] = Min_N_imm;

    // continuous inputs
    if (i < drv.inputRates.size()) {
      const auto& in = drv.inputRates[i];
      for (int k = 0; k < NUM_POOLS; k++) d[k] += in[k];
      fr[i][F_INPUT_C] = totalC(in);
      fr[i][F_INPUT_N] = totalN(in);
    }
  }

  // signed transport of dissolved LMW-C and LMW-N across the lower boundary of each layer
  for (size_t i = 0; i < n; i++) {
    const double flux = i < drv.waterFlux_mm.size() ? drv.waterFlux_mm[i] : 0.0;
    const double wUp = drv.layers[i].water_kg_ha;
    double tC, tN;
    if (i + 1 < n) {
      const double wLow = drv.layers[i + 1].water_kg_ha;
      tC = soluteTransport(flux, y[i][LMW_C], wUp, &y[i + 1][LMW_C], wLow, 0, _p.leachCap);
      tN = soluteTransport(flux, y[i][LMW_N], wUp, &y[i + 1][LMW_N], wLow, 0, _p.leachCap);
      dy[i + 1][LMW_C] += tC;
      dy[i + 1][LMW_N] += tN;
    } else {
      tC = soluteTransport(flux, y[i][LMW_C], wUp, nullptr, 0, drv.deepConcC, _p.leachCap);
      tN = soluteTransport(flux, y[i][LMW_N], wUp, nullptr, 0, drv.deepConcN, _p.leachCap);
    }
    dy[i][LMW_C] -= tC;
    dy[i][LMW_N] -= tN;
    fr[i][F_LEACH_C_OUT] = tC;
    fr[i][F_LEACH_N_OUT] = tN;
  }
}

void Profile::stepRK4(std::vector<LayerPools>& y, const std::function<const Drivers&(double)>& drivers,
                      std::vector<LayerFluxes>& fluxes, double t0, double dt, bool* negativeStage) const {
  const size_t n = y.size();
  std::vector<LayerPools> k1, k2, k3, k4, tmp(n);
  std::vector<LayerFluxes> f1, f2, f3, f4;
  bool neg = false;
  auto stage = [&](const std::vector<LayerPools>& k, double w) {
    for (size_t i = 0; i < n; i++) {
      for (int j = 0; j < NUM_POOLS; j++) {
        tmp[i][j] = y[i][j] + w * dt * k[i][j];
        if (tmp[i][j] < 0 && y[i][j] >= 0) neg = true;
      }
    }
  };

  derivatives(y, drivers(t0), k1, f1);
  stage(k1, 0.5);
  derivatives(tmp, drivers(t0 + 0.5 * dt), k2, f2);
  stage(k2, 0.5);
  derivatives(tmp, drivers(t0 + 0.5 * dt), k3, f3);
  stage(k3, 1.0);
  derivatives(tmp, drivers(t0 + dt), k4, f4);

  if (fluxes.size() != n) fluxes.assign(n, LayerFluxes{});
  for (size_t i = 0; i < n; i++) {
    for (int k = 0; k < NUM_POOLS; k++) {
      const double before = y[i][k];
      y[i][k] = y[i][k] + dt * (k1[i][k] + 2 * k2[i][k] + 2 * k3[i][k] + k4[i][k]) / 6.0;
      if (y[i][k] < 0 && before >= 0) neg = true;
    }
    for (int k = 0; k < NUM_FLUXES; k++) {
      fluxes[i][k] += dt * (f1[i][k] + 2 * f2[i][k] + 2 * f3[i][k] + f4[i][k]) / 6.0;
    }
  }
  if (negativeStage) *negativeStage = *negativeStage || neg;
}

void Profile::stepRK4(std::vector<LayerPools>& y, const Drivers& d, std::vector<LayerFluxes>& fluxes) const {
  stepRK4(y, [&d](double) -> const Drivers& { return d; }, fluxes);
}

int Profile::stepAdaptive(std::vector<LayerPools>& y, const std::function<const Drivers&(double)>& drivers,
                          std::vector<LayerFluxes>& fluxes) const {
  const size_t nl = y.size();
  if (fluxes.size() != nl) fluxes.assign(nl, LayerFluxes{});
  const int maxN = std::max(1, int(_p.maxSubsteps));
  if (_p.maxTurnoverPerStep <= 0) {
    stepRK4(y, drivers, fluxes);
    return 1;
  }
  // largest relative net loss of a pool per day at the start of the day
  std::vector<LayerPools> dy;
  std::vector<LayerFluxes> fr;
  derivatives(y, drivers(0.0), dy, fr);
  double turnover = 0;
  for (size_t i = 0; i < nl; i++) {
    for (int k = 0; k < NUM_POOLS; k++) {
      if (y[i][k] > 1e-6 && dy[i][k] < 0) turnover = std::max(turnover, -dy[i][k] / y[i][k]);
    }
  }
  int n = std::min(maxN, std::max(1, int(std::ceil(turnover / _p.maxTurnoverPerStep - 1e-12))));
  while (true) {
    auto yt = y;
    std::vector<LayerFluxes> ft(nl, LayerFluxes{});
    bool neg = false;
    const double dt = 1.0 / n;
    for (int s = 0; s < n; s++) stepRK4(yt, drivers, ft, s * dt, dt, &neg);
    if (neg && n < maxN) {
      n = std::min(maxN, 2 * n);
      continue;
    }
    y = yt;
    for (size_t i = 0; i < nl; i++) for (int k = 0; k < NUM_FLUXES; k++) fluxes[i][k] += ft[i][k];
    return n;
  }
}

int Profile::stepAdaptive(std::vector<LayerPools>& y, const Drivers& d, std::vector<LayerFluxes>& fluxes) const {
  return stepAdaptive(y, [&d](double) -> const Drivers& { return d; }, fluxes);
}

std::pair<double, double> Profile::balanceResidual(const std::vector<LayerPools>& before,
                                                   const std::vector<LayerPools>& after,
                                                   const std::vector<LayerFluxes>& f) {
  double dC = 0, dN = 0;
  for (size_t i = 0; i < before.size(); i++) {
    dC += totalC(after[i]) - totalC(before[i]) + f[i][F_CO2] - f[i][F_INPUT_C];
    dN += totalN(after[i]) - totalN(before[i]) + f[i][F_MIN_GROSS] - f[i][F_MIN_IMM] - f[i][F_INPUT_N];
  }
  if (!before.empty()) {
    dC += f.back()[F_LEACH_C_OUT];
    dN += f.back()[F_LEACH_N_OUT];
  }
  return {dC, dN};
}

} // namespace samm
