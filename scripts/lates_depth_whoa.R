# examining read depth for lates filtered VCFs
library(tidyverse)

depth <- read_table2("../results/lates_all_092320_rawvariants_depth_byInd2.txt",
                     col_names=TRUE)
depth2 <- depth %>%
  rowwise() %>%
  mutate(mean_dp = mean(CEW16_135_2:JAR18_140026, na.rm=T)) %>%
  dplyr::select(c(chrom,pos,overall_dp,mean_dp))

depth %>%
  filter(overall_dp > 5 & overall_dp < 10000) %>%
  ggplot(aes(), alpha=0.7) +
    geom_histogram(aes(x=overall_dp),breaks=c(0,2,5,10,20,30,40,50,seq(100,2000,by=50)),closed="left")


## using whoa package
## analyzing all species combined + each species separately
## using miss0.5 maf0.01 thin90 datasets
library(whoa)
source("theme_custom.R")
library(patchwork)

lates_vcf <- vcfR::read.vcfR("lates_all_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
gfreqs <- exp_and_obs_geno_freqs(lates_vcf)
p.lates <- geno_freqs_scatter(gfreqs,max_plot_loci=100000) +
  theme_custom()
binned_lates <- infer_m(lates_vcf, minBin=20000)

lsta_vcf <- vcfR::read.vcfR("lsta_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
gfreqs_lsta <- exp_and_obs_geno_freqs(lsta_vcf)
p.lsta <- geno_freqs_scatter(gfreqs_lsta,max_plot_loci=100000) +
  theme_custom()
binned_lsta <- infer_m(lsta_vcf, minBin=2000)

lmic_vcf <- vcfR::read.vcfR("lmic_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
gfreqs_lmic <- exp_and_obs_geno_freqs(lmic_vcf)
p.lmic <- geno_freqs_scatter(gfreqs_lmic,max_plot_loci=100000) +
  theme_custom()
binned_lmic <- infer_m(lmic_vcf, minBin=2000)

lmar_vcf <- vcfR::read.vcfR("lmar_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
gfreqs_lmar <- exp_and_obs_geno_freqs(lmar_vcf)
p.lmar <- geno_freqs_scatter(gfreqs_lmar,max_plot_loci=100000) +
  theme_custom()
binned_lmar <- infer_m(lmar_vcf, minBin=2000)

lang_vcf <- vcfR::read.vcfR("lang_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
gfreqs_lang <- exp_and_obs_geno_freqs(lang_vcf)
p.lang <- geno_freqs_scatter(gfreqs_lang,max_plot_loci=100000) +
  theme_custom()
binned_lang <- infer_m(lang_vcf, minBin=2000)

p2.lates <- posteriors_plot(binned_lates$m_posteriors) +
  theme_custom() + theme(legend.position="none") + 
  ylab("Posterior mean estimate\nof miscall rate")
p2.lsta <- posteriors_plot(binned_lsta$m_posteriors) + 
  theme_custom() + theme(legend.position="none") + 
  ylab("Posterior mean estimate\nof miscall rate")
p2.lmic <- posteriors_plot(binned_lmic$m_posteriors) + 
  theme_custom() + theme(legend.position="none") + 
  ylab("Posterior mean estimate\nof miscall rate")
p2.lmar <- posteriors_plot(binned_lmar$m_posteriors) + 
  theme_custom() + theme(legend.position="none") + 
  ylab("Posterior mean estimate\nof miscall rate")
p2.lang <- posteriors_plot(binned_lang$m_posteriors) + 
  theme_custom() + theme(legend.position="none") + 
  ylab("Posterior mean estimate\nof miscall rate")

p.lates + p.lsta + p.lmic + p.lmar + p.lang +
  p2.lates + p2.lsta + p2.lmic + p2.lmar + p2.lang +
  plot_layout(widths=c(2,1),ncol=2,byrow=FALSE)
