# Litter chemistry for SAMM

`litter-chemistry.csv` holds lignin and polyphenol contents (kg per kg dry matter) of all crop residues
(`monica-parameters/crop-residues`, looked up by `species`) and organic fertilisers
(`monica-parameters/organic-fertilisers`, looked up by `id`), with separate values for dead roots of each crop.
MONICA's SAMM module converts them to the structural (STR) share of the added carbon:

    STR share of C = (0.63 * lignin + 0.53 * polyphenols) / C content

After editing the CSV, regenerate the C++ table and rebuild MONICA:

    python3 data/samm/generate_litter_chemistry.py

Values in a residue or fertiliser JSON file override the table: `SAMM_Lignin`, `SAMM_Polyphenols`,
`SAMM_RootLignin`, `SAMM_RootPolyphenols` (kg kg-1 DM), and for organic fertilisers optionally
`SAMM_ShareToLMW`, `SAMM_ShareToMAO` (share of the added C and N going directly to these pools).

## How the values were derived

Every row says in `basis_shoot` / `basis_root` where its value comes from:

- **measured**: the source, plant part, harvest stage and number of observations (n) are given;
- **estimate**: derived from a related material or the crop group, as described in the row. These should be
  replaced by measured values where available (for example from the experiment being simulated).

Lignin is Van Soest acid detergent lignin throughout (as in the SAMM paper). Polyphenol (soluble, Folin)
contents of shoots and manures are group estimates; only the roots of 20 crops have measured values (R14).
Polyphenols have a small weight in the structural share (0.5-4% of DM against 3-25% lignin).

| Code | Source | Used for |
|---|---|---|
| T21 | Thiébeau, P., Jensen, L.S., Ferchaud, F., Recous, S. (2021): Dataset of biomass and chemical quality of crop residues from European areas. Data in Brief 37, 107227, https://doi.org/10.1016/j.dib.2021.107227. Data: https://doi.org/10.15454/LBI3U7 (version 3, mean data files) | shoot lignin of field crops, cover crops, forage crops and vegetables; mature roots of maize, rape, ryegrass and brassica vegetables |
| R14 | Redin, M., Guénon, R., Recous, S., Schmatz, R., de Freitas, L.L., Aita, C., Giacomini, S.J. (2014): Carbon mineralization in soil of roots from twenty crop species, as affected by their chemical composition and botanical family. Plant and Soil 378, 205-214, https://doi.org/10.1007/s11104-013-2021-5 (Table 2) | root lignin and polyphenols of 20 crops, family means for the other crops |
| E24 | EOM4SOIL - Physico-chemical characteristics of external organic matters (EOMs) database, version 1 (2024-10-25), https://doi.org/10.5281/zenodo.13969793; Michaud, A.M. et al.: Dataset on physico-chemical characteristics of Exogenous Organic Matters (EOMs) gathered from various European countries, Data in Brief | lignin of manures, slurries, composts and sewage sludge (French data; lignin in % of organic matter multiplied by the organic matter content of the same EOM type) |

Caveats:

- R14 sampled roots at flowering; the authors note that root lignin at maturity is likely higher
  (Abiven et al. 2005 report 160-260 g/kg). Where T21 has mature roots (maize 7.7%, rape 15.7%) those are used.
- T21 cover crops are immature (autumn destruction); MONICA's cover crop residue files are used for such material.
- Manure lignin varies widely with bedding and storage; deep-litter manures and slurries without own data were
  assigned the value of the closest manure type.
- The Khon Kaen materials of the SAMM paper (rice straw, groundnut, tamarind, dipterocarp) are not in MONICA's
  parameter files; their values are in `tests/samm/khonkaen_monica/make_setup.py` (Laub et al. 2024, Table 1).
