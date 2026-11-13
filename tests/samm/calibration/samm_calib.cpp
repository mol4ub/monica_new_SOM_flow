// Calibration of SAMM for the Khon Kaen experiment with a real-valued genetic algorithm, using the objective
// function of the SAMM paper (Laub et al. 2024, GMD 17:931, section 2.4): the Gaussian log-likelihood of litterbag C,
// microbial N, SOC, aggregate C and free MAOC, without covariances, with the median measurement variance of each
// data type. The model is the SAMM C++ core (src/core/samm.cpp) with the set-up of Samm_main.R (one 0-15 cm layer,
// HYDRUS drivers, litter additions, daily root inputs, litterbag pools), integrated with deSolve's classical RK4.
// Litterbag: the bag added on day 4414 is followed past the next addition (see bagDeriv), unlike Samm_main.R v1.1.
//
// build: g++ -O3 -std=c++17 -pthread -I src/core tests/samm/calibration/samm_calib.cpp src/core/samm.cpp -o samm_calib
// usage: samm_calib <climate_file.csv> <AllMeasurements.csv> <mode> [options]
//   modes:  eval <params.csv> <out.csv>    run one parameter set, write daily pools of all treatments
//           ga <out_prefix>                genetic algorithm, writes <prefix>_best.csv and <prefix>_log.csv
//           sir <prior.csv> <samples.csv>  draws --n parameter sets from the prior (name,mean,sd,lo,hi: normal, truncated
//                                          to [lo,hi]) and writes each set with its log-likelihood (resampling:
//                                          sir_resample.py)
//           dream <prior.csv> <start.csv> <chains.csv>  DREAM(ZS) MCMC (Vrugt 2016) with --chains chains (default 4)
//                                          started around <start.csv>; target = prior x likelihood; writes every
//                                          --thin-th generation of all chains (summary: dream_summary.py)
//           band <posterior.csv> <out.csv> runs --n parameter sets drawn from a posterior sample file (parameter
//                                          columns by name) and writes daily 2.5/50/97.5% quantiles per treatment
//   options: --modifiers SAMM|KirschbaumMONICA (default KirschbaumMONICA), --threads n, --pop n, --gens n,
//            --seed n, --thetaSat x (default 0.468), --moistOpt x (default 0.45), --n n (default 200000),
//            --chains n (4), --maxGen n (60000), --thin n (10), --rhat x (stop when all R-hat < x, default 1.1),
//            --lik paper|ar1 (default paper; see Config::ar1), --exclude "<data type>" (ar1 only, repeatable)
#include "samm.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace samm;

struct Day { double temp, theta, leach_mm, pF; };
struct Obs { int trt, day, type; double value; };
enum ObsType { LITTERBAG, MICN, SOC, AGGC, MAOC, MICC, SOILN, NTYPES };
static const char* typeNames[NTYPES] = {"Litterbag C", "Microbial N", "SOC", "Aggregate C", "Free MOAC", "Microbial C",
                                        "Soil N"};
// short names for the error-model hyperparameters sigma_<key> and tau_<key> (likelihood "ar1")
static const char* typeKeys[NTYPES] = {"Litterbag", "MicN", "SOC", "AggC", "MAOC", "MicC", "SoilN"};
static const int NPM = 23;  // number of model parameters incl. observation model (the hyperparameters follow them)
static const double FN_USED = 3.1;  // factor used to convert ninhydrin-reactive N to microbial N in the data
static int hyperIndex[NTYPES] = {-1, -1, -1, -1, -1, -1, -1};  // position of sigmaRatio_<type> in the parameters

static std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> f; std::stringstream ss(s); std::string x;
  while (std::getline(ss, x, sep)) f.push_back(x);
  return f;
}

static std::vector<Day> readClimate(const std::string& path) {
  std::ifstream in(path); std::string line; std::getline(in, line);
  std::vector<Day> days;
  while (std::getline(in, line)) {
    auto f = split(line, ';');
    if (f.size() < 6 || f[2].empty() || f[5].empty()) break;
    days.push_back({std::stod(f[2]), std::stod(f[3]), std::stod(f[4]), std::stod(f[5])});
  }
  return days;
}

static const char* trtNames[5] = {"CT", "DP", "GN", "RS", "TM"};
struct Treatment { double LitterC, LitterN, LitterStructuralPercent; };
static const Treatment trts[5] = {{0, 0, 0}, {4530, 57, 0.32}, {3880, 228, 0.13}, {3670, 47, 0.06}, {4300, 136, 0.17}};
static const std::set<int> addDays = {1, 382, 759, 1125, 1490, 1874, 2231, 2596, 2961, 3323, 3688, 4049, 4414,
                                      4779, 5144, 5509, 5874, 6240, 6605, 6970, 7335, 7702, 8068, 8432, 8797, 9162};
static const int LAST = 9070;
static const int LITTERBAG_DAY = 4414;  // litter addition of the litterbag experiment

// calibrated parameters (Table 4 of the paper) and their search bounds; upper bounds of k_MAO (10 x SAMM v1.1),
// CUE_STR, MicMaxCN and Daily_Litter_C (weed roots) limit the GA to plausible values (without them it ran to
// k_MAO 0.01, CUE_STR 0.99, MicMaxCN 15 and 8 kg C/ha/d of root input)
struct PDef { const char* name; double lo, hi; };
static std::vector<PDef> pdefs = {
  {"k_STR", 0.0015, 0.08}, {"k_LAB", 0.02, 0.9}, {"k_MIC", 0.003, 0.2}, {"k_MAO", 0.0002, 0.0041},
  {"mu_max", 0.05, 0.95}, {"kAgg", 0.015, 0.8}, {"K_M_MIC_C", 5, 400}, {"m_MIC", 0.0, 0.02},
  {"K_LMW_MAO", 0.005, 0.5}, {"CUE_STR", 0.1, 0.7}, {"CUE_LAB", 0.1, 1.0}, {"MicMinCN", 3.0, 7.5},
  {"MicMaxCN", 7.5, 12.0}, {"f_MICMAOM", 0.02, 0.8}, {"pc_STR_LAB", 0.5, 8.0}, {"aggfactSTR_C", 0.0, 5.0},
  {"aggfactMAO_C", 0.0, 8.0}, {"NonMicAgg", 0.0, 150.0}, {"Daily_Litter_C", 0.5, 5.0},
  {"Daily_Litter_CN", 20.0, 300.0}, {"Daily_LitterSTR_C", 0.01, 0.6}, {"betaMIC", 1.0, 2.0},
  // observation model: true factor microbial N / ninhydrin-reactive N; the data were converted with 3.1
  // (Amato & Ladd 1988; Laub et al. 2022), so reported microbial N = model MIC_N * 3.1 / fN_ninhydrin
  {"fN_ninhydrin", 2.0, 5.0}};
static int NP = int(pdefs.size());

struct Config {
  bool sammModifiers = false;
  //! false ("paper"): independent Gaussian errors with the median measurement variance per data type, data types
  //! litterbag C, microbial N, SOC, aggregate C, free MAOC (Laub et al. 2024).
  //! true ("ar1"): all seven data types (except --exclude); per treatment and data type, residual = independent
  //! measurement error (median measurement variance var_type) + model error with standard deviation
  //! sigmaRatio_<type> * sqrt(var_type) and exponential temporal correlation exp(-dt / tau_<type>); exact likelihood
  //! with a Kalman filter.
  bool ar1 = false;
  //! data types left out of the likelihood "ar1" (--exclude "<type name>")
  std::set<int> excluded;
  double thetaSat = 0.468, moistOpt = 0.45;
};

static Params toParams(const std::vector<double>& x) {
  Params p = Params::samm_v11();
  p.withoutMineralN();
  p.maxTurnoverPerStep = 0;  // one RK4 step per day, as deSolve::rk4
  double* m[] = {&p.k_STR, &p.k_LAB, &p.k_MIC, &p.k_MAO, &p.mu_max, &p.kAgg, &p.K_M_MIC_C, &p.m_MIC, &p.K_LMW_MAO,
                 &p.CUE_STR, &p.CUE_LAB, &p.MicMinCN, &p.MicMaxCN, &p.f_MICMAOM, &p.pc_STR_LAB, &p.aggfactSTR_C,
                 &p.aggfactMAO_C, &p.NonMicAgg};
  for (int i = 0; i < 18; i++) *m[i] = x[i];
  p.betaMIC = x[21];
  return p;
}

// litterbag pools of the R code: decompose like STR and LAB, are not protected in aggregates, emptied 13 days
// before each addition and refilled with the added litter
struct Bag { double STR, LAB, LABN; };

struct Result { std::vector<std::array<double, NUM_POOLS + 3>> rows; bool ok = true; };

// runs one treatment; rows[t-1] = state at the start of day t (R output row "time = t")
static Result run(const std::vector<double>& x, int ti, const std::vector<Day>& clim, const Config& cfg, bool keepRows,
                  const std::vector<int>* needDays, std::vector<std::array<double, NUM_POOLS + 3>>* atDays) {
  const Params p = toParams(x);
  const double DLC = x[18], DLCN = x[19], DLSTR = x[20];
  LayerEnv env; env.thickness = 0.15; env.bulkDensity = 1450; env.siltClay = 0.176;
  Profile prof(p, {env});
  const auto& t = trts[ti];

  LayerPools y0{};
  y0[LMW_C] = 100; y0[MIC_C] = 100; y0[MAO_C] = 1665 - 575;
  y0[AggSTR_C] = 1400 * 0.1; y0[AggLAB_C] = 1400 * 0.9; y0[AggMAO_C] = 575;
  y0[LMW_N] = 5; y0[MIC_N] = 10;
  const double initLitter = 4900 - y0[LMW_C] - y0[MIC_C] - y0[MAO_C] - y0[AggMAO_C] - y0[AggSTR_C] - y0[AggLAB_C];
  y0[STR_C] = initLitter * DLSTR; y0[LAB_C] = initLitter - y0[STR_C]; y0[LAB_N] = y0[LAB_C] / DLCN;
  y0[MAO_N] = y0[MAO_C] / 9; y0[AggLAB_N] = y0[AggLAB_C] / DLCN; y0[AggMAO_N] = y0[AggMAO_C] / 9;
  Bag b0{0.1, 0.1, 0.01};

  std::vector<LayerPools> y = {y0};
  Bag bag = b0;
  Drivers d; d.layers.resize(1); d.waterFlux_mm.resize(1); d.inputRates.resize(1);
  auto setDrivers = [&](double T) -> int {
    const int idx = int(std::floor(T));
    const auto& c = clim.at(idx - 1);
    double cf;
    if (cfg.sammModifiers) cf = rateModWaterSAMM(c.pF) * rateModTempSAMM(c.temp);
    else cf = rateModTempKirschbaum(c.temp) * std::exp(-18.0 * std::pow(c.theta / cfg.thetaSat - cfg.moistOpt, 2));
    d.layers[0].corrF = cf;
    d.layers[0].water_kg_ha = (c.theta * 0.15 * 10000) * 1000;
    d.waterFlux_mm[0] = c.leach_mm;
    auto& in = d.inputRates[0]; in.fill(0);
    if (addDays.count(idx)) { in[STR_C] += t.LitterC * t.LitterStructuralPercent;
      in[LAB_C] += t.LitterC * (1 - t.LitterStructuralPercent); in[LAB_N] += t.LitterN; }
    if (idx > 1) { in[STR_C] += DLC * DLSTR; in[LAB_C] += DLC * (1 - DLSTR); in[LAB_N] += DLC / DLCN; }
    return idx;
  };
  auto bagDeriv = [&](const Bag& bg, double micC, int idx, double cf) -> Bag {
    const double aMIC = std::max(micC / (p.K_M_MIC_C + micC), p.a_MIC_min);
    double prot = bg.STR * p.pc_STR_LAB / bg.LAB; if (prot > 1) prot = 1;
    Bag db{-bg.STR * p.k_STR * aMIC * cf, -bg.LAB * (1 - prot) * p.k_LAB * aMIC * cf,
           -bg.LABN * (1 - prot) * p.k_LAB * aMIC * cf};
    // the bag of the litterbag experiment (added on day LITTERBAG_DAY) is followed until the last measurement
    // (day 4799, after the next addition): from then on the bag is neither emptied nor refilled (Samm_main.R v1.1
    // empties and refills every year; the paper version had second-year litterbag pools for this)
    const bool bagFrozen = idx > LITTERBAG_DAY;
    if (!bagFrozen && addDays.count(idx + 13)) { db = {-bg.STR * 0.99999, -bg.LAB * 0.99999, -bg.LABN * 0.99999}; }
    if (!bagFrozen && addDays.count(idx)) { db.STR += t.LitterC * t.LitterStructuralPercent;
      db.LAB += t.LitterC * (1 - t.LitterStructuralPercent); db.LABN += t.LitterN; }
    return db;
  };

  Result res;
  auto record = [&](int time) {
    std::array<double, NUM_POOLS + 3> r;
    for (int k = 0; k < NUM_POOLS; k++) r[k] = y[0][k];
    r[NUM_POOLS] = bag.STR; r[NUM_POOLS + 1] = bag.LAB; r[NUM_POOLS + 2] = bag.LABN;
    if (keepRows) res.rows.push_back(r);
    if (atDays) for (size_t i = 0; i < needDays->size(); i++) if ((*needDays)[i] == time) (*atDays)[i] = r;
  };
  record(1);
  std::vector<LayerPools> k1, k2, k3, k4, tmp(1);
  std::vector<LayerFluxes> f;
  for (int time = 1; time < LAST; time++) {
    Bag kb[4];
    const double offs[4] = {0, 0.5, 0.5, 1};
    const double w[4] = {0.5, 0.5, 1, 0};
    std::vector<LayerPools>* ks[4] = {&k1, &k2, &k3, &k4};
    tmp = y; Bag bt = bag;
    for (int s = 0; s < 4; s++) {
      const int idx = setDrivers(time + offs[s]);
      prof.derivatives(tmp, d, *ks[s], f);
      kb[s] = bagDeriv(bt, tmp[0][MIC_C], idx, d.layers[0].corrF);
      if (s < 3) {
        for (int k = 0; k < NUM_POOLS; k++) tmp[0][k] = y[0][k] + w[s] * (*ks[s])[0][k];
        bt = {bag.STR + w[s] * kb[s].STR, bag.LAB + w[s] * kb[s].LAB, bag.LABN + w[s] * kb[s].LABN};
      }
    }
    for (int k = 0; k < NUM_POOLS; k++) {
      y[0][k] = y[0][k] + (k1[0][k] + 2 * k2[0][k] + 2 * k3[0][k] + k4[0][k]) / 6.0;
      if (!std::isfinite(y[0][k]) || y[0][k] < -1e-6) res.ok = false;
    }
    bag.STR += (kb[0].STR + 2 * kb[1].STR + 2 * kb[2].STR + kb[3].STR) / 6.0;
    bag.LAB += (kb[0].LAB + 2 * kb[1].LAB + 2 * kb[2].LAB + kb[3].LAB) / 6.0;
    bag.LABN += (kb[0].LABN + 2 * kb[1].LABN + 2 * kb[2].LABN + kb[3].LABN) / 6.0;
    if (!res.ok && !keepRows) return res;
    record(time + 1);
  }
  return res;
}

static double simValue(const std::array<double, NUM_POOLS + 3>& r, int type) {
  switch (type) {
    case LITTERBAG: return r[NUM_POOLS] + r[NUM_POOLS + 1];
    case MICN: return r[MIC_N];
    // SOC: all pools except the free (non-aggregated) litter pools STR_C and LAB_C (paper, Appendix A2)
    case SOC: return r[LMW_C] + r[MIC_C] + r[MAO_C] + r[AggSTR_C] + r[AggLAB_C] + r[AggMAO_C];
    case AGGC: return r[AggSTR_C] + r[AggLAB_C] + r[AggMAO_C];
    case MAOC: return r[MAO_C];
    case MICC: return r[MIC_C];
    // soil N: the N of the same pools as SOC (without free LAB_N)
    case SOILN: return r[LMW_N] + r[MIC_N] + r[MAO_N] + r[AggLAB_N] + r[AggMAO_N];
  }
  return 0;
}

struct Data {
  std::vector<Obs> obs;
  double var[NTYPES];
  std::vector<int> days[5];  // observation days per treatment
  std::vector<std::vector<size_t>> series;  // observation indices per treatment and data type, sorted by day
};

static Data readData(const std::string& path) {
  Data D;
  std::ifstream in(path); std::string line; std::getline(in, line);
  std::vector<double> vars[NTYPES];
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto f = split(line, ';');
    if (f.size() < 8) continue;
    int type = -1;
    for (int k = 0; k < NTYPES; k++) if (f[3] == typeNames[k]) type = k;
    int trt = -1;
    for (int k = 0; k < 5; k++) if (f[2] == trtNames[k]) trt = k;
    if (type < 0 || trt < 0) continue;
    D.obs.push_back({trt, std::stoi(f[0]), type, std::stod(f[6])});
    vars[type].push_back(std::stod(f[7]));
  }
  for (int k = 0; k < NTYPES; k++) {
    auto v = vars[k]; std::sort(v.begin(), v.end());
    size_t n = v.size();
    D.var[k] = n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  }
  for (const auto& o : D.obs) D.days[o.trt].push_back(std::max(o.day, 1));
  for (int t = 0; t < 5; t++) for (int ty = 0; ty < NTYPES; ty++) {
    std::vector<size_t> idx;
    for (size_t i = 0; i < D.obs.size(); i++) if (D.obs[i].trt == t && D.obs[i].type == ty) idx.push_back(i);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return D.obs[a].day < D.obs[b].day; });
    if (!idx.empty()) D.series.push_back(idx);
  }
  return D;
}

// objective: Gaussian log-likelihood (constant terms included, as the paper's AIC uses the likelihood)
static double logLik(const std::vector<double>& x, const std::vector<Day>& clim, const Data& D, const Config& cfg,
                     std::vector<double>* sims = nullptr) {
  double ll = 0;
  std::vector<double> sim(D.obs.size(), std::numeric_limits<double>::quiet_NaN());
  for (int ti = 0; ti < 5; ti++) {
    std::vector<std::array<double, NUM_POOLS + 3>> at(D.days[ti].size());
    auto r = run(x, ti, clim, cfg, false, &D.days[ti], &at);
    if (!r.ok) return -std::numeric_limits<double>::infinity();
    size_t j = 0;
    for (size_t i = 0; i < D.obs.size(); i++) if (D.obs[i].trt == ti) {
      sim[i] = simValue(at[j++], D.obs[i].type);
      if (D.obs[i].type == MICN) sim[i] *= FN_USED / x[22];  // microbial N as reported with the factor 3.1
    }
  }
  if (sims) *sims = sim;
  if (!cfg.ar1) {
    for (size_t i = 0; i < D.obs.size(); i++) {
      const auto& o = D.obs[i];
      if (o.type == MICC || o.type == SOILN) continue;  // not used in the paper's calibration
      ll += -0.5 * std::log(2 * M_PI * D.var[o.type]) - 0.5 * (o.value - sim[i]) * (o.value - sim[i]) / D.var[o.type];
    }
    return ll;
  }
  // Kalman filter per treatment and data type: latent model error m_t (stationary, variance sigma^2, correlation
  // exp(-dt/tau)) plus independent measurement error with the median measurement variance of the data type
  for (const auto& series : D.series) {
    const int ty = D.obs[series[0]].type;
    if (hyperIndex[ty] < 0) continue;  // excluded data type
    const double ratio = x[hyperIndex[ty]], tau = x[hyperIndex[ty] + 1], r = D.var[ty], q = ratio * ratio * r;
    double m = 0, P = q; int prevDay = 0;
    for (size_t k = 0; k < series.size(); k++) {
      const auto& o = D.obs[series[k]];
      if (k > 0) { const double phi = std::exp(-(o.day - prevDay) / tau); m *= phi; P = phi * phi * P + q * (1 - phi * phi); }
      const double e = o.value - sim[series[k]], S = P + r, v = e - m;
      ll += -0.5 * (std::log(2 * M_PI * S) + v * v / S);
      const double K = P / S; m += K * v; P *= (1 - K);
      prevDay = o.day;
    }
  }
  return ll;
}

static std::vector<double> readParams(const std::string& path) {
  std::map<std::string, double> m;
  std::ifstream in(path); std::string line;
  while (std::getline(in, line)) { auto f = split(line, ','); if (f.size() >= 2 && f[0] != "name") m[f[0]] = std::stod(f[1]); }
  std::vector<double> x(NP);
  for (int i = 0; i < NP; i++) {
    if (!m.count(pdefs[i].name) && std::string(pdefs[i].name) == "betaMIC") m["betaMIC"] = 1.0;  // older files: linear death
    if (!m.count(pdefs[i].name) && std::string(pdefs[i].name) == "fN_ninhydrin") m["fN_ninhydrin"] = FN_USED;  // no correction
    if (!m.count(pdefs[i].name)) throw std::runtime_error(std::string("missing parameter ") + pdefs[i].name);
    x[i] = m[pdefs[i].name];
  }
  return x;
}

static void writeParams(const std::string& path, const std::vector<double>& x, double ll) {
  std::ofstream out(path); out.precision(10);
  out << "name,value\n";
  for (int i = 0; i < NP; i++) out << pdefs[i].name << "," << x[i] << "\n";
  out << "logLik," << ll << "\n";
}

int main(int argc, char** argv) {
  if (argc < 4) { std::cerr << "usage: samm_calib climate.csv AllMeasurements.csv eval|ga ...\n"; return 1; }
  const auto clim = readClimate(argv[1]);
  const auto D = readData(argv[2]);
  const std::string mode = argv[3];
  Config cfg; int threads = std::max(1u, std::thread::hardware_concurrency()), pop = 100, gens = 300;
  unsigned seed = 1;
  long nSir = 200000;
  int nChains = 4, thin = 10; long maxGen = 60000; double rhatStop = 1.1;
  std::vector<std::string> pos;
  for (int a = 4; a < argc; a++) {
    std::string s = argv[a];
    if (s == "--modifiers") cfg.sammModifiers = std::string(argv[++a]) == "SAMM";
    else if (s == "--threads") threads = std::stoi(argv[++a]);
    else if (s == "--pop") pop = std::stoi(argv[++a]);
    else if (s == "--gens") gens = std::stoi(argv[++a]);
    else if (s == "--seed") seed = std::stoul(argv[++a]);
    else if (s == "--thetaSat") cfg.thetaSat = std::stod(argv[++a]);
    else if (s == "--moistOpt") cfg.moistOpt = std::stod(argv[++a]);
    else if (s == "--n") nSir = std::stol(argv[++a]);
    else if (s == "--chains") nChains = std::stoi(argv[++a]);
    else if (s == "--maxGen") maxGen = std::stol(argv[++a]);
    else if (s == "--thin") thin = std::stoi(argv[++a]);
    else if (s == "--rhat") rhatStop = std::stod(argv[++a]);
    else if (s == "--lik") cfg.ar1 = std::string(argv[++a]) == "ar1";
    else if (s == "--exclude") { std::string t = argv[++a]; bool found = false;
      for (int k = 0; k < NTYPES; k++) if (t == typeNames[k]) { cfg.excluded.insert(k); found = true; }
      if (!found) throw std::runtime_error("unknown data type for --exclude: " + t); }
    else pos.push_back(s);
  }
  std::fprintf(stderr, "median variances:");
  for (int k = 0; k < NTYPES; k++) std::fprintf(stderr, " %s %.1f", typeNames[k], D.var[k]);
  std::fprintf(stderr, " (%zu observations)\n", D.obs.size());
  static std::vector<std::string> hyperNames;
  if (cfg.ar1) {
    hyperNames.reserve(2 * NTYPES);
    for (int k = 0; k < NTYPES; k++) {
      if (cfg.excluded.count(k)) continue;
      hyperNames.push_back(std::string("sigmaRatio_") + typeKeys[k]); hyperNames.push_back(std::string("tau_") + typeKeys[k]);
      hyperIndex[k] = int(pdefs.size());
      pdefs.push_back({hyperNames[hyperNames.size() - 2].c_str(), 0.0, 3.0});
      pdefs.push_back({hyperNames[hyperNames.size() - 1].c_str(), 1.0, 365.0});
    }
    NP = int(pdefs.size());
    std::fprintf(stderr, "likelihood ar1: %d data series, %d parameters incl. error model\n", int(D.series.size()), NP);
  }

  if (mode == "eval") {
    auto x = readParams(pos.at(0));
    std::vector<double> sims;
    double ll = logLik(x, clim, D, cfg, &sims);
    std::printf("logLik %.3f\n", ll);
    std::ofstream out(pos.at(1)); out.precision(12);
    out << "time";
    for (int k = 0; k < NUM_POOLS; k++) out << "," << poolName(k);
    out << ",LitterBagSTR_C,LitterBagLAB_C,LitterBagLAB_N,TRT\n";
    for (int ti = 0; ti < 5; ti++) {
      auto r = run(x, ti, clim, cfg, true, nullptr, nullptr);
      for (size_t i = 0; i < r.rows.size(); i++) {
        out << i + 1; for (double v : r.rows[i]) out << "," << v; out << "," << trtNames[ti] << "\n";
      }
    }
    if (pos.size() > 2) {
      std::ofstream o2(pos[2]); o2.precision(10);
      o2 << "TRT,day,type,obs,sim\n";
      for (size_t i = 0; i < D.obs.size(); i++)
        o2 << trtNames[D.obs[i].trt] << "," << D.obs[i].day << "," << typeNames[D.obs[i].type] << "," << D.obs[i].value << "," << sims[i] << "\n";
    }
    return 0;
  }

  if (mode == "ga") {
    // real-valued GA in the spirit of R's GA package: linear-rank selection, blend crossover, random mutation, elitism
    const std::string prefix = pos.at(0);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    std::vector<std::vector<double>> P(pop, std::vector<double>(NP));
    for (auto& ind : P) for (int i = 0; i < NP; i++) ind[i] = pdefs[i].lo + U(rng) * (pdefs[i].hi - pdefs[i].lo);
    for (size_t s = 0; s + 1 < pos.size() && int(s) < pop; s++) {  // starting points, clipped into the bounds
      P[s] = readParams(pos[s + 1]);
      for (int i = 0; i < NP; i++) P[s][i] = std::min(pdefs[i].hi, std::max(pdefs[i].lo, P[s][i]));
    }
    std::vector<double> fit(pop);
    auto evalAll = [&](std::vector<std::vector<double>>& Q, std::vector<double>& F, int from) {
      std::atomic<int> next(from);
      std::vector<std::thread> th;
      for (int k = 0; k < threads; k++) th.emplace_back([&]() {
        for (int i; (i = next++) < int(Q.size());) F[i] = logLik(Q[i], clim, D, cfg);
      });
      for (auto& t : th) t.join();
    };
    evalAll(P, fit, 0);
    std::ofstream log(prefix + "_log.csv");
    log << "generation,best,median\n";
    const int nElite = std::max(2, pop / 20);
    std::vector<double> best; double bestF = -INFINITY; int stall = 0;
    for (int g = 0; g < gens; g++) {
      std::vector<int> ord(pop); for (int i = 0; i < pop; i++) ord[i] = i;
      std::sort(ord.begin(), ord.end(), [&](int a, int b) { return fit[a] > fit[b]; });
      if (fit[ord[0]] > bestF + 1e-6) { bestF = fit[ord[0]]; best = P[ord[0]]; stall = 0; } else stall++;
      std::vector<double> fs; for (double v : fit) if (std::isfinite(v)) fs.push_back(v);
      std::sort(fs.begin(), fs.end());
      log << g << "," << bestF << "," << (fs.empty() ? NAN : fs[fs.size() / 2]) << "\n"; log.flush();
      std::fprintf(stderr, "gen %3d best %.3f median %.3f\n", g, bestF, fs.empty() ? NAN : fs[fs.size() / 2]);
      if (g % 10 == 0) writeParams(prefix + "_best.csv", best, bestF);
      if (stall >= 60) break;
      // linear rank selection probabilities (best gets 2/pop, worst 0)
      std::vector<double> cum(pop); double c = 0;
      for (int r = 0; r < pop; r++) { c += 2.0 * (pop - r) / (pop * (pop + 1.0)); cum[r] = c; }
      auto pick = [&]() { double u = U(rng) * c; int r = int(std::lower_bound(cum.begin(), cum.end(), u) - cum.begin());
                          return P[ord[std::min(r, pop - 1)]]; };
      std::vector<std::vector<double>> Q;
      for (int e = 0; e < nElite; e++) Q.push_back(P[ord[e]]);
      const double mutScale = std::max(0.02, 0.25 * (1.0 - double(g) / gens));
      while (int(Q.size()) < pop) {
        auto a = pick(), b = pick();
        std::vector<double> ch(NP);
        for (int i = 0; i < NP; i++) {
          double lo = std::min(a[i], b[i]), hi = std::max(a[i], b[i]), dd = hi - lo;
          double v = (U(rng) < 0.8) ? lo - 0.5 * dd + U(rng) * 2 * dd : a[i];  // BLX-0.5 crossover
          if (U(rng) < 0.1) {                                                    // mutation
            std::normal_distribution<double> N(0, mutScale * (pdefs[i].hi - pdefs[i].lo));
            v += N(rng);
          }
          ch[i] = std::min(pdefs[i].hi, std::max(pdefs[i].lo, v));
        }
        Q.push_back(ch);
      }
      std::vector<double> F(pop);
      for (int e = 0; e < nElite; e++) F[e] = fit[ord[e]];
      evalAll(Q, F, nElite);
      P = Q; fit = F;
    }
    writeParams(prefix + "_best.csv", best, bestF);
    std::printf("best logLik %.3f\n", bestF);
    return 0;
  }

  if (mode == "sir") {
    // prior: normal distributions truncated to [lo, hi] (rejection sampling); MicMinCN < MicMaxCN
    std::map<std::string, std::array<double, 4>> pr;
    { std::ifstream in(pos.at(0)); std::string line;
      while (std::getline(in, line)) { auto f = split(line, ','); if (f.size() >= 5 && f[0] != "name")
        pr[f[0]] = {std::stod(f[1]), std::stod(f[2]), std::stod(f[3]), std::stod(f[4])}; } }
    for (int i = 0; i < NP; i++) if (!pr.count(pdefs[i].name)) throw std::runtime_error(std::string("prior missing for ") + pdefs[i].name);
    std::mt19937_64 rng(seed);
    std::vector<std::vector<double>> X(nSir, std::vector<double>(NP));
    for (long n = 0; n < nSir; n++) {
      do {
        for (int i = 0; i < NP; i++) {
          const auto& q = pr[pdefs[i].name];
          std::normal_distribution<double> N(q[0], q[1]);
          double v; do { v = N(rng); } while (v < q[2] || v > q[3]);
          X[n][i] = v;
        }
      } while (X[n][11] >= X[n][12]);  // MicMinCN < MicMaxCN
    }
    std::vector<double> L(nSir);
    std::atomic<long> next(0), done(0);
    std::vector<std::thread> th;
    for (int k = 0; k < threads; k++) th.emplace_back([&]() {
      for (long i; (i = next++) < nSir;) {
        L[i] = logLik(X[i], clim, D, cfg);
        long d = ++done;
        if (d % 10000 == 0) std::fprintf(stderr, "%ld of %ld simulations\n", d, nSir);
      }
    });
    for (auto& t : th) t.join();
    std::ofstream out(pos.at(1)); out.precision(10);
    for (int i = 0; i < NP; i++) out << pdefs[i].name << ",";
    out << "logLik\n";
    for (long n = 0; n < nSir; n++) { for (int i = 0; i < NP; i++) out << X[n][i] << ","; out << L[n] << "\n"; }
    return 0;
  }

  if (mode == "dream") {
    // ---- prior (normal, truncated to [lo, hi]) ----
    std::vector<std::array<double, 4>> pr(NP);
    { std::map<std::string, std::array<double, 4>> m; std::ifstream in(pos.at(0)); std::string line;
      while (std::getline(in, line)) { auto f = split(line, ','); if (f.size() >= 5 && f[0] != "name")
        m[f[0]] = {std::stod(f[1]), std::stod(f[2]), std::stod(f[3]), std::stod(f[4])}; }
      for (int i = 0; i < NP; i++) { if (!m.count(pdefs[i].name)) throw std::runtime_error(std::string("prior missing for ") + pdefs[i].name);
                                     pr[i] = m[pdefs[i].name]; } }
    auto inBounds = [&](const std::vector<double>& x) {
      for (int i = 0; i < NP; i++) if (x[i] < pr[i][2] || x[i] > pr[i][3]) return false;
      return x[11] < x[12];
    };
    auto logPrior = [&](const std::vector<double>& x) -> double {
      if (!inBounds(x)) return -INFINITY;
      double lp = 0; for (int i = 0; i < NP; i++) { double z = (x[i] - pr[i][0]) / pr[i][1]; lp += -0.5 * z * z; }
      return lp;
    };
    auto fold = [&](std::vector<double>& x) {  // fold back into [lo, hi] (DREAM boundary handling)
      for (int i = 0; i < NP; i++) {
        const double lo = pr[i][2], hi = pr[i][3], w = hi - lo;
        for (int k = 0; k < 10 && (x[i] < lo || x[i] > hi); k++) {
          if (x[i] < lo) x[i] = hi - (lo - x[i]); else x[i] = lo + (x[i] - hi);
        }
        if (x[i] < lo || x[i] > hi) x[i] = lo + 0.5 * w;
      }
    };
    const auto x0 = readParams(pos.at(1));
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    std::normal_distribution<double> N01(0, 1);
    const int d = NP, nc = nChains, delta = 1, kArch = 10;
    // ---- initial archive Z: 10 d points around the start (sd = 0.1 x prior sd), chains start from the first nc ----
    std::vector<std::vector<double>> Z;
    while (int(Z.size()) < 10 * d) {
      std::vector<double> x = x0;
      if (!Z.empty()) for (int i = 0; i < d; i++) x[i] += 0.1 * pr[i][1] * N01(rng);
      fold(x); if (inBounds(x)) Z.push_back(x);
    }
    std::vector<std::vector<double>> X(Z.begin(), Z.begin() + nc);
    std::vector<double> LL(nc), LP(nc);
    auto evalPar = [&](std::vector<std::vector<double>>& xs, std::vector<double>& ll) {
      std::vector<std::thread> th;
      for (size_t c = 0; c < xs.size(); c++) th.emplace_back([&, c]() { ll[c] = logLik(xs[c], clim, D, cfg); });
      for (auto& t : th) t.join();
    };
    evalPar(X, LL);
    for (int c = 0; c < nc; c++) LP[c] = LL[c] + logPrior(X[c]);
    std::ofstream out(pos.at(2)); out.precision(10);
    out << "gen,chain"; for (int i = 0; i < NP; i++) out << "," << pdefs[i].name; out << ",logLik,logPost\n";
    std::vector<std::vector<std::vector<double>>> hist(nc);  // thinned chain states for R-hat
    const double CRs[3] = {1.0 / 3, 2.0 / 3, 1.0};
    long accepted = 0, proposed = 0;
    for (long g = 1; g <= maxGen; g++) {
      std::vector<std::vector<double>> P(nc); std::vector<double> jac(nc, 0.0);
      for (int c = 0; c < nc; c++) {
        std::vector<double> x = X[c], dx(d, 0.0);
        if (U(rng) < 0.1) {
          // snooker update
          const auto& z = Z[size_t(U(rng) * Z.size())];
          const auto& z1 = Z[size_t(U(rng) * Z.size())]; const auto& z2 = Z[size_t(U(rng) * Z.size())];
          double nn = 0; std::vector<double> dir(d); for (int i = 0; i < d; i++) { dir[i] = x[i] - z[i]; nn += dir[i] * dir[i]; }
          if (nn > 0) {
            double p1 = 0, p2 = 0; for (int i = 0; i < d; i++) { p1 += z1[i] * dir[i]; p2 += z2[i] * dir[i]; }
            const double gs = 1.2 + U(rng);
            for (int i = 0; i < d; i++) dx[i] = gs * (p1 - p2) / nn * dir[i];
            std::vector<double> xn(d); double nn2 = 0;
            for (int i = 0; i < d; i++) { xn[i] = x[i] + dx[i]; nn2 += (xn[i] - z[i]) * (xn[i] - z[i]); }
            jac[c] = 0.5 * (d - 1) * (std::log(std::max(nn2, 1e-300)) - std::log(nn));
          }
        } else {
          // parallel direction update with randomised subspace sampling
          const double CR = CRs[int(U(rng) * 3) % 3];
          std::vector<int> A; for (int i = 0; i < d; i++) if (U(rng) < CR) A.push_back(i);
          if (A.empty()) A.push_back(int(U(rng) * d) % d);
          const double gamma = (U(rng) < 0.2) ? 1.0 : 2.38 / std::sqrt(2.0 * delta * A.size());
          for (int k = 0; k < delta; k++) {
            const auto& z1 = Z[size_t(U(rng) * Z.size())]; const auto& z2 = Z[size_t(U(rng) * Z.size())];
            for (int i : A) dx[i] += gamma * (1 + 0.05 * (2 * U(rng) - 1)) * (z1[i] - z2[i]);
          }
          for (int i : A) dx[i] += 1e-6 * std::abs(x[i] + 1e-12) * N01(rng);
        }
        for (int i = 0; i < d; i++) x[i] += dx[i];
        fold(x); P[c] = x;
      }
      std::vector<double> PL(nc);
      evalPar(P, PL);
      for (int c = 0; c < nc; c++) {
        const double lp = PL[c] + logPrior(P[c]);
        proposed++;
        if (std::isfinite(lp) && std::log(U(rng)) < lp - LP[c] + jac[c]) { X[c] = P[c]; LL[c] = PL[c]; LP[c] = lp; accepted++; }
      }
      if (g % kArch == 0) for (int c = 0; c < nc; c++) Z.push_back(X[c]);
      if (g % thin == 0) for (int c = 0; c < nc; c++) {
        out << g << "," << c; for (double v : X[c]) out << "," << v; out << "," << LL[c] << "," << LP[c] << "\n";
        hist[c].push_back(X[c]);
      }
      if (g % 1000 == 0) {
        out.flush();
        // Gelman-Rubin R-hat on the second half of the thinned chains
        const size_t n = hist[0].size(), h = n / 2, m = n - h;
        double rmax = 0; int worst = 0;
        for (int i = 0; i < d && m > 2; i++) {
          std::vector<double> mean(nc), var(nc); double gm = 0;
          for (int c = 0; c < nc; c++) { double s = 0; for (size_t t = h; t < n; t++) s += hist[c][t][i]; mean[c] = s / m; gm += mean[c] / nc;
                                          double v = 0; for (size_t t = h; t < n; t++) v += (hist[c][t][i] - mean[c]) * (hist[c][t][i] - mean[c]); var[c] = v / (m - 1); }
          double B = 0, W = 0; for (int c = 0; c < nc; c++) { B += (mean[c] - gm) * (mean[c] - gm); W += var[c]; }
          B *= double(m) / (nc - 1); W /= nc;
          const double r = W > 0 ? std::sqrt(((m - 1.0) / m * W + B / m) / W) : 1e9;
          if (r > rmax) { rmax = r; worst = i; }
        }
        const double bestLL = *std::max_element(LL.begin(), LL.end());
        std::fprintf(stderr, "gen %6ld  acceptance %.3f  max R-hat %.3f (%s)  logLik of chains %.1f %.1f %.1f %.1f\n",
                     g, double(accepted) / proposed, rmax, pdefs[worst].name, LL[0], nc > 1 ? LL[1] : NAN, nc > 2 ? LL[2] : NAN,
                     nc > 3 ? LL[3] : NAN);
        (void)bestLL;
        if (g >= 10000 && rmax < rhatStop) { std::fprintf(stderr, "converged (all R-hat < %.2f)\n", rhatStop); break; }
      }
    }
    return 0;
  }

  if (mode == "band") {
    std::ifstream in(pos.at(0)); std::string line; std::getline(in, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    auto hdr = split(line, ','); std::vector<int> col(NP, -1);
    for (int i = 0; i < NP; i++) for (size_t k = 0; k < hdr.size(); k++) if (hdr[k] == pdefs[i].name) col[i] = int(k);
    std::vector<std::vector<double>> all;
    while (std::getline(in, line)) {
      auto f = split(line, ',');
      std::vector<double> x(NP);
      for (int i = 0; i < NP; i++) x[i] = std::stod(f.at(col[i]));
      all.push_back(x);
    }
    std::mt19937_64 rng(seed); std::vector<std::vector<double>> S;
    const long n = std::min<long>(nSir, long(all.size()));
    for (long k = 0; k < n; k++) S.push_back(all[size_t(std::uniform_real_distribution<double>(0, 1)(rng) * all.size())]);
    const char* qn[7] = {"SOC", "Microbial N", "Microbial C", "Aggregate C", "Free MOAC", "Litterbag C", "Soil N"};
    auto qv = [](const std::array<double, NUM_POOLS + 3>& r, int q) {
      switch (q) { case 0: return r[LMW_C] + r[MIC_C] + r[MAO_C] + r[AggSTR_C] + r[AggLAB_C] + r[AggMAO_C]; case 1: return r[MIC_N];
        case 2: return r[MIC_C]; case 3: return r[AggSTR_C] + r[AggLAB_C] + r[AggMAO_C]; case 4: return r[MAO_C];
        case 5: return r[NUM_POOLS] + r[NUM_POOLS + 1];
        default: return r[LMW_N] + r[MIC_N] + r[MAO_N] + r[AggLAB_N] + r[AggMAO_N]; } };
    std::ofstream out(pos.at(1)); out.precision(8);
    out << "TRT,time,quantity,q025,q50,q975\n";
    for (int ti = 0; ti < 5; ti++) {
      // vals[q][day][k]
      std::vector<std::vector<std::vector<float>>> vals(7, std::vector<std::vector<float>>(LAST, std::vector<float>(n)));
      std::atomic<long> next(0); std::vector<std::thread> th;
      for (int t = 0; t < threads; t++) th.emplace_back([&]() {
        for (long k; (k = next++) < n;) {
          auto r = run(S[k], ti, clim, cfg, true, nullptr, nullptr);
          for (int day = 0; day < LAST && day < int(r.rows.size()); day++) for (int q = 0; q < 7; q++) vals[q][day][k] = float(qv(r.rows[day], q));
        }
      });
      for (auto& t : th) t.join();
      for (int day = 0; day < LAST; day++) for (int q = 0; q < 7; q++) {
        auto v = vals[q][day]; std::sort(v.begin(), v.end());
        auto Q = [&](double p) { return v[std::min<size_t>(v.size() - 1, size_t(p * v.size()))]; };
        out << trtNames[ti] << "," << day + 1 << "," << qn[q] << "," << Q(0.025) << "," << Q(0.5) << "," << Q(0.975) << "\n";
      }
    }
    return 0;
  }
  std::cerr << "unknown mode " << mode << "\n";
  return 1;
}
