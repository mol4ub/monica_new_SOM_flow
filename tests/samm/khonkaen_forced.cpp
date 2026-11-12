// Runs the SAMM C++ core (src/core/samm.cpp) for the Khon Kaen litter experiment with exactly the
// set-up of the SAMM v1.1 R code (tests/samm/data/Samm_main.R): one 0-15 cm layer, soil drivers from
// climate_file.csv, litter added as in the R code and integrated with deSolve's classical RK4.
//
// build: g++ -O2 -std=c++17 -I src/core tests/samm/khonkaen_forced.cpp src/core/samm.cpp -o khonkaen_forced
// run:   ./khonkaen_forced tests/samm/data/climate_file.csv out.csv [mineralN]
//        (mineralN: also switch on the scenario D extensions, f_min/f_turn/k_des, with no mineral N pool)

#include "samm.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace samm;

struct Day { double temp, theta, leach_mm, pF; };

static std::vector<Day> readClimate(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(in, line); // header
  std::vector<Day> days;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string f;
    std::vector<std::string> fields;
    while (std::getline(ss, f, ';')) fields.push_back(f);
    if (fields.size() < 6 || fields[2].empty() || fields[5].empty()) break; // trailing empty rows
    days.push_back({std::stod(fields[2]), std::stod(fields[3]), std::stod(fields[4]), std::stod(fields[5])});
  }
  return days;
}

struct Treatment { const char* name; double LitterC, LitterN, LitterStructuralPercent; };

int main(int argc, char** argv) {
  if (argc < 3) { std::cerr << "usage: khonkaen_forced climate_file.csv out.csv [mineralN]\n"; return 1; }
  const bool mineralN = argc > 3 && std::string(argv[3]) == "mineralN";
  const auto climate = readClimate(argv[1]);

  Params p = Params::samm_v11();
  if (!mineralN) p.withoutMineralN();
  const double simulated_depth = 0.15, bulk_density = 1.45, si_cl = 0.176;
  LayerEnv env;
  env.thickness = simulated_depth;
  env.bulkDensity = bulk_density * 1000;
  env.siltClay = si_cl;
  Profile profile(p, {env});

  // daily root inputs of the R code (calibrated, unavoidable plant growth)
  const double Daily_Litter_C = 2.28, Daily_Litter_CN = 142.9, Daily_LitterSTR_C = 0.16;
  const std::set<int> addDays = {1, 382, 759, 1125, 1490, 1874, 2231, 2596, 2961, 3323, 3688, 4049, 4414,
                                 4779, 5144, 5509, 5874, 6240, 6605, 6970, 7335, 7702, 8068, 8432, 8797, 9162};
  const Treatment trts[] = {{"CT", 0, 0, 0}, {"DP", 4530, 57, 0.32}, {"GN", 3880, 228, 0.13},
                            {"RS", 3670, 47, 0.06}, {"TM", 4300, 136, 0.17}};

  std::ofstream out(argv[2]);
  out << "time";
  for (int k = 0; k < NUM_POOLS; k++) out << "," << poolName(k);
  out << ",CO2_cum,LeachC_cum,LeachN_cum,MinGross_cum,TRT\n";
  out.precision(17);

  for (const auto& t : trts) {
    // initial state as in the R code
    LayerPools y0{};
    const double InitialSystemC = 4900;
    y0[LMW_C] = 100; y0[MIC_C] = 100; y0[MAO_C] = 1665 - 575;
    y0[AggSTR_C] = 1400 * 0.1; y0[AggLAB_C] = 1400 * 0.9; y0[AggMAO_C] = 575;
    y0[LMW_N] = 5; y0[MIC_N] = 10;
    double InitialLitterC = InitialSystemC - y0[LMW_C] - y0[MIC_C] - y0[MAO_C] - y0[AggMAO_C] - y0[AggSTR_C] - y0[AggLAB_C];
    y0[STR_C] = InitialLitterC * Daily_LitterSTR_C;
    y0[LAB_C] = InitialLitterC - y0[STR_C];
    y0[LAB_N] = y0[LAB_C] / Daily_Litter_CN;
    y0[MAO_N] = y0[MAO_C] / 9;
    y0[AggLAB_N] = y0[AggLAB_C] / Daily_Litter_CN;
    y0[AggMAO_N] = y0[AggMAO_C] / 9;

    std::vector<LayerPools> y = {y0};
    std::vector<LayerFluxes> fluxes(1, LayerFluxes{});
    Drivers d;
    d.layers.resize(1);
    d.waterFlux_mm.resize(1);
    d.inputRates.resize(1);

    // drivers at time T, indexed by floor(T) as Weather[floor(Time)] in R (1-based)
    auto drivers = [&](double T) -> const Drivers& {
      int idx = int(std::floor(T));
      const auto& c = climate.at(idx - 1);
      d.layers[0].corrF = rateModWaterSAMM(c.pF) * rateModTempSAMM(c.temp);
      d.layers[0].water_kg_ha = (c.theta * simulated_depth * 10000) * 1000;
      d.waterFlux_mm[0] = c.leach_mm;
      auto& in = d.inputRates[0];
      in.fill(0);
      if (addDays.count(idx)) {
        in[STR_C] += t.LitterC * t.LitterStructuralPercent;
        in[LAB_C] += t.LitterC * (1 - t.LitterStructuralPercent);
        in[LAB_N] += t.LitterN;
      }
      if (idx > 1) {
        in[STR_C] += Daily_Litter_C * Daily_LitterSTR_C;
        in[LAB_C] += Daily_Litter_C * (1 - Daily_LitterSTR_C);
        in[LAB_N] += Daily_Litter_C / Daily_Litter_CN;
      }
      return d;
    };

    double maxResC = 0, maxResN = 0;
    auto write = [&](int time) {
      out << time;
      for (int k = 0; k < NUM_POOLS; k++) out << "," << y[0][k];
      out << "," << fluxes[0][F_CO2] << "," << fluxes[0][F_LEACH_C_OUT] << "," << fluxes[0][F_LEACH_N_OUT]
          << "," << fluxes[0][F_MIN_GROSS] << "," << t.name << "\n";
    };
    write(1);
    for (int time = 1; time < 9070; time++) {
      auto before = y;
      std::vector<LayerFluxes> daily(1, LayerFluxes{});
      profile.stepRK4(y, [&](double off) -> const Drivers& { return drivers(time + off); }, daily);
      auto res = Profile::balanceResidual(before, y, daily);
      maxResC = std::max(maxResC, std::abs(res.first));
      maxResN = std::max(maxResN, std::abs(res.second));
      for (int k = 0; k < NUM_FLUXES; k++) fluxes[0][k] += daily[0][k];
      write(time + 1);
    }
    std::printf("%s: max daily balance residual C %.2e N %.2e kg/ha\n", t.name, maxResC, maxResN);
  }
  return 0;
}
