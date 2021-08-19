# lates-popgen
[ Note: repository is a work in progress, June 2021 ] Repository for scripts associated with Rick, Junker, et al.'s manuscript on population genetics in Lake Tanganyika's Lates spp.

## Contents
### Data
- lates\_all\_metadata.csv : metadata for individuals included in the GBS dataset
- lates\_fst\_by\_sampling\_site.csv : Fst estimates by species between all sampling sites
- lates\_metadata\_assignments.csv : metadata including entropy cluster assignments for the GBS dataset
- lates\_rad\_metadata.csv : metadata for individuals included in the RAD dataset
- lktang\_samplingLocs.csv : latitude and longitudes for named sampling sites
- Note: filtered VCFs can be found in the associated Zenodo repository at http://doi.org/10.5281/zenodo.5216259
### Scripts
*Scripts for command line genomic analyses*
- checkHetsIndvVCF.sh
- run_angsd_thetas.sh
- run_angsd_thetas_array.sh
- Entropy scripts:
  - vcf2mpgl.R
  - plot_entropy_results.R
*Scripts for analysis in R*
- latesIBD_lang.R
- latesIBD_lmar.R
- latesIBD_lmic.R
- latesIBD_lsta.R
- lates_depth_whoa.R
- lates_fst_ind.R
- plot_emu_pca.R
*Helper scripts*
- packages_funcs.R
- theme_custom.R
