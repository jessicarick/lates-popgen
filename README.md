# lates-popgen
[ Note: repository is a work in progress, June 2021 ] Repository for scripts associated with Rick, Junker, et al.'s manuscript on population genetics in Lake Tanganyika's Lates spp.

## Contents
### Data
- lates\_all\_metadata.csv : metadata for individuals included in the GBS dataset
- lates\_fst\_by\_sampling\_site.csv : Fst estimates by species between all sampling sites
- lates\_metadata\_assignments.csv : metadata including entropy cluster assignments for the GBS dataset
- lates\_rad\_metadata.csv : metadata for individuals included in the RAD dataset
- lktang\_samplingLocs.csv : latitude and longitudes for named sampling sites
- Note: filtered VCFs can be found in the associated Zenodo repository at http://doi.org/10.5281/zenodo.5216259 and fastq files are uploaded to NCBI SRA as a part of Bioproject PRJNA776855 (will be released upon publication of the manuscript)
### Scripts
- bash/cluster:
  - checkHetsIndvVCF.sh : script for checking ratio of reads supporting minor vs major allele at heterozygous sites
  - run_angsd_thetas.sh : script for running ANGSD to calculate genetic diversity within species
  - run_angsd_thetas_array.sh : script for running ANGSD to calculate genetic diversity within species, parallelized using SLURM arrays
- entropy-related scripts
  - vcf2mpgl.R : conversion script modified from [the entropy vignette](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/auxfiles/inputdataformat.R)
  - run_entropy_array.sh : script to run entropy as a bash array
  - run_entropy_parallel_array.sh : script to run entropy chains in parallel and with different k values in a bash array
  - lates_entropy_plots.R
  - lates_entropy_plots_perSpp.R
- R scripts:
  - latesIBD_\*.R : scripts for calculating site-wise Fsts and testing for isolation-by-distance in each species
  - lates_depth_whoa.R : script for assessing heterozygote excess as it relates to read depth
  - lates_fst_ind.R : script for computing individual-based Fsts and analyzing patterns of isolation-by-distance, including testing hypotheses related to IBD patterns
  - lates_fst_pop.R : script for calculating FST between species and between populations within species, as well as heterozygosity within populations
  - packages_funcs.R : helper script containing packages and functions needed for R analyses
  - plot_emu_pca.R : script for plotting results from EMU PCA analyses
  - theme_custom.R : helper script specifying base ggplot theme
  - vcf2mpgl.R : script for converting VCF to mpgl for input into entropy
