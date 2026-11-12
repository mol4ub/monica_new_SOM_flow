/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
Built-in litter chemistry (lignin and polyphenol contents) of MONICA's crop residues and organic
fertilisers for the SAMM soil organic matter module. The values and their sources are in
data/samm/litter-chemistry.csv; samm-litter-chemistry.cpp is generated from it by
data/samm/generate_litter_chemistry.py.
*/

#pragma once

#include <string>

namespace samm {

struct LitterChemistry {
  const char* key;          //!< normalised species name (crop residues) or id (organic fertilisers)
  double lignin;            //!< [kg kg-1 DM]
  double polyphenols;       //!< [kg kg-1 DM]
  double rootLignin;        //!< [kg kg-1 DM]
  double rootPolyphenols;   //!< [kg kg-1 DM]
};

//! lower case, only letters and digits ("Sudan gras" -> "sudangras", "bell - sweet pepper" -> "bellsweetpepper")
std::string normaliseChemistryKey(const std::string& s);

//! nullptr if the species / id is not in the table
const LitterChemistry* residueChemistry(const std::string& species);
const LitterChemistry* fertiliserChemistry(const std::string& id);

//! used when a material is not in the table
const LitterChemistry& defaultResidueChemistry();
const LitterChemistry& defaultFertiliserChemistry();

//! carbon contents converting lignin and polyphenol mass to structural C (as in the SAMM paper's set-up)
constexpr double C_CONTENT_LIGNIN = 0.63;
constexpr double C_CONTENT_POLYPHENOLS = 0.53;

} // namespace samm
