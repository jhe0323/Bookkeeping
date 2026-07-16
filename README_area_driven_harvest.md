# Area-driven LUH2 harvest patch

This patch changes wood harvest to a pure LUH2-area-driven scheme.

Changed files:
- `events.py`
- `LULCCSimulator.py`

Main behavior:
1. Only LUH2 `*_harv` variables are used.
2. LUH2 `*_bioh` variables are ignored and are not read by `_parse_wood_harvest()`.
3. Harvested biomass is computed from the model state on the LUH2 harvest area:
   `removed biomass = current modeled biomass × LUH2 harvest area share`.
4. No expansion outside LUH2 harvest footprint is allowed.
   Therefore `harvest_forced_biomass`, `harvest_extra_area`, `harvest_unmet_biomass`, and `harvest_unmet_raw_biomass` remain zero in normal area-driven runs.
5. `harvest_requested_biomass` is kept for output compatibility, but now means the area-implied modeled biomass removed, not LUH2 `*_bioh` demand.

To use it, replace the corresponding files in your code directory with the two files in this folder.
