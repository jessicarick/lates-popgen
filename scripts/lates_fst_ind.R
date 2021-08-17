###################
## Script for ind-level fst analyses
##
## Analysis for Rick et al., J Hered
## Written by J. Rick, jrick@uwyo.edu
## Last update: Summer 2021
###################

##################### 
## Loading packages and data
#####################
## Load required packages
source("packages_funcs.R")
library(broom)

vcf.rename <- function(x) {
  col.names <- unlist(strsplit(indNames(x),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
  col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
  col.names.clean[col.names.clean == "CEW16_135_2"] <- "CEW16_135"
  indNames(x) <- col.names.clean
  ploidy(x) <- 2
  x <- gl.compliance.check(x)
  return(x)
}

## Importing the SNP data and metadata
## working with the datasets created with the following filters (from VCFtools): 
## missing data allowed = 50%, MAF = 0.01 cutoff, min depth = 5, thinned 90bp

lates_vcfR <- read.vcfR("../../lates_popgen_github_copy/data/lates_all_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lang_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lang_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmar_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lmar_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmic_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lmic_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lsta_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_rad_vcfR <- read.vcfR("../data/combined_lsta_noHets3-4_variants_0.5_maf0.01_dp5_thin90.recode.vcf")

fishinfo <- read.csv('../data/lates_all_metadata.csv',
                     header=TRUE,stringsAsFactors = FALSE) %>%
  mutate(Moran_FishID = ind,
         SL_mm = as.numeric(SL_mm),
         TL_mm = as.numeric(TL_mm))

locs <- read.csv("../data/lktang_samplingLocs.csv",header=TRUE,row.names = 1)
colnames(locs) <- c("Site","Latitude","Longitude","coords")
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

Dgeo <- geodist::geodist(locs, paired = FALSE, sequential = FALSE, pad = FALSE,
                              measure = "geodesic")
row.names(Dgeo) <- locs$Site
colnames(Dgeo) <- locs$Site
Dgeo_long <- subset(melt(Dgeo), value!=0)

#------------------------------#
## LOADING VCFS, FIXING NAMES, AND CALCULATING INDIVIDUAL PAIRWISE FSTS
#------------------------------#

lang_gen <- vcf.rename(vcfR2genlight(lang_vcfR))
pop(lang_gen) <- indNames(lang_gen)
lang_fst_ind <- reich.fst(lang_gen,bootstrap=FALSE,plot=FALSE,verbose=TRUE)

lmar_gen <- vcf.rename(vcfR2genlight(lmar_vcfR))
pop(lmar_gen) <- indNames(lmar_gen)
lmar_fst_ind <- reich.fst(lmar_gen,bootstrap=FALSE,plot=FALSE,verbose=TRUE)

lmic_gen <- vcf.rename(vcfR2genlight(lmic_vcfR))
pop(lmic_gen) <- indNames(lmic_gen)
lmic_fst_ind  <- reich.fst(lmic_gen,bootstrap=FALSE,plot=FALSE,verbose=TRUE)

lsta_gen <- vcf.rename(vcfR2genlight(lsta_vcfR))
pop(lsta_gen) <- indNames(lsta_gen)
lsta_fst_ind  <- reich.fst(lsta_gen,bootstrap=FALSE,plot=FALSE,verbose=TRUE)

#------------------------------#
## ANALYZING INDIVIDUAL-LEVEL FST
#------------------------------#

## combining data into one dataframe
lang_fst_ind_2 <- lang_fst_ind$fsts %>%
  as_data_frame() %>%
  mutate(pop1 = colnames(.)) %>%
  pivot_longer(cols=!starts_with("pop"),names_to="pop2",values_to="fst") %>%
  mutate(spp = "lang") %>%
  drop_na()
lmar_fst_ind_2 <- lmar_fst_ind$fsts %>%
  as_data_frame() %>%
  mutate(pop1 = colnames(.)) %>%
  pivot_longer(cols=!starts_with("pop"),names_to="pop2",values_to="fst") %>%
  mutate(spp = "lmar") %>%
  drop_na()
lmic_fst_ind_2 <- lmic_fst_ind$fsts %>%
  as_data_frame() %>%
  mutate(pop1 = colnames(.)) %>%
  pivot_longer(cols=!starts_with("pop"),names_to="pop2",values_to="fst") %>%
  mutate(spp = "lmic") %>%
  drop_na()
lsta_fst_ind_2 <- lsta_fst_ind$fsts %>%
  as_data_frame() %>%
  mutate(pop1 = colnames(.)) %>%
  pivot_longer(cols=!starts_with("pop"),names_to="pop2",values_to="fst") %>%
  mutate(spp = "lsta") %>%
  drop_na()

lates_fst_ind_all <- lang_fst_ind_2 %>%
  bind_rows(lmar_fst_ind_2,lmic_fst_ind_2,lsta_fst_ind_2) %>%
  left_join(fishinfo,by=c("pop1" = "ind")) %>%
  left_join(fishinfo,by=c("pop2" = "ind"),suffix=c(".pop1",".pop2")) %>%
  mutate(diff_sl = abs(as.numeric(SL_mm.pop1)-as.numeric(SL_mm.pop2)),
         juv = case_when(juvenile.pop1 == "Y" & juvenile.pop2 == "Y" ~ "juv-juv",
                         juvenile.pop1 == "Y" & juvenile.pop2 == "N" ~ "juv-adult",
                         juvenile.pop1 == "N" & juvenile.pop2 == "Y" ~ "juv-adult",
                         juvenile.pop1 == "N" & juvenile.pop2 == "N" ~ "adult-adult",
                         juvenile.pop1 == "" | juvenile.pop2 == "" ~ ""),
         sampling_season.pop1 = case_when(sampling_month.pop1 == "August" ~ "Fall",
                                          sampling_month.pop1 == "September" ~ "Fall",
                                          sampling_month.pop1 == "May" ~ "Spring",
                                          sampling_month.pop1 == "July" ~ "Summer",
                                          sampling_month.pop1 == "" ~ ""),
         sampling_season.pop2 = case_when(sampling_month.pop2 == "August" ~ "Fall",
                                          sampling_month.pop2 == "September" ~ "Fall",
                                          sampling_month.pop2 == "May" ~ "Spring",
                                          sampling_month.pop2 == "July" ~ "Summer",
                                          sampling_month.pop2 == "" ~ ""),
         combined_month = paste0(sampling_month.pop1,"-",sampling_month.pop2),
         combined_season = case_when(sampling_season.pop1 == sampling_season.pop2 ~ sampling_season.pop1,
                                     sampling_season.pop1 != sampling_season.pop2 ~ paste0(sampling_season.pop1,"-",sampling_season.pop2)),
         combined_season = case_when(combined_season %in% c("Fall-Spring","Spring-Fall") ~ "Spring-Fall",
                                     combined_season %in% c("Fall-Summer", "Summer-Fall") ~ "Summer-Fall",
                                     combined_season %in% c("Spring-Summer","Summer-Spring") ~ "Spring-Summer",
                                     TRUE ~ combined_season)) %>%
  left_join(Dgeo_long,by=c("sampling_loc.pop1" = "Var1","sampling_loc.pop2" = "Var2")) %>%
  left_join(Dgeo_long,by=c("sampling_loc.pop1" = "Var2", "sampling_loc.pop2" = "Var1")) %>%
  mutate(geo_dist_km = case_when(value.x > 0 ~ value.x,
                                 value.y > 0 ~ value.y)) 

## Q1. Do juveniles show more IBD than adults?
## Checking and plotting IBD slopes by juv categories
lates_fst_ind_all %>%
  group_by(spp,juv) %>%
  do(fitJuv = tidy(lm((fst/(1-fst)) ~ log(geo_dist_km), data=.))) %>%
  unnest(fitJuv)

lates_fst_ind_all %>% 
  filter(!is.na(fst) & !is.na(geo_dist_km)) %>%
  ggplot() +
    geom_point(aes(x=log(geo_dist_km),y=(fst/(1-fst)),col=juv,fill=juv)) +
    geom_smooth(aes(x=log(geo_dist_km),y=(fst/(1-fst)),col=juv),method="lm",lty=2,lwd=2) +
    theme_custom() +
    scale_fill_npg() +
    facet_wrap(~spp, scales="free")


## Q2. Are fish of similar sizes more closely related to one another?
lates_fst_ind_all %>%
  group_by(spp,juv) %>%
  do(fitSL = tidy(lm(fst ~ diff_sl, data=.))) %>%
  unnest(fitSL)

lates_fst_ind_all %>% 
  filter(!is.na(fst) 
           #sampling_month.pop1 != "" & 
           #sampling_month.pop2 != "" & 
           #juv == "juv-juv"
           #sampling_year.pop1 %in% c("2017","2018") &
           #sampling_year.pop2 %in% c("2017","2018")
  ) %>%
  ggplot() +
    geom_point(aes(x=diff_sl,y=fst)) +
    geom_smooth(aes(x=diff_sl,y=fst),method="lm",lty=2,lwd=2,col="black") +
    theme_custom() +
    scale_fill_npg() +
    facet_wrap(~spp, scales="free")


## Q3. Do fish of the same juvenile status or same size belong to the same entropy clusters?
## NOTE: NOT COMPLETED YET
lates_fst_ind_all %>%
  left_join(entropy_assignments) %>% # dataframe from plot_entropy_results.R script
  filter(!is.na(fst)) %>%
  ggplot() +
    geom_jitter(aes(x=k2,y=diff_sl,fill=juv),size=4,pch=21, alpha=0.7,width=0.1,height=0) +
    geom_boxplot(aes(x=k2,y=diff_sl),fill=NA,outlier.color=NA) +
    theme_custom() +
    scale_fill_npg() +
    facet_wrap(~spp,scales="free")

## Q4. Are FST estimates for juveniles higher than those for adult fish?
lang_dat <- lates_fst_ind_all %>%
  filter(!is.na(fst) & sampling_loc.pop1 == sampling_loc.pop2) %>%
  filter(spp == "lang")

t.test(lang_dat$fst[lang_dat$juv == "juv-juv"],lang_dat$fst[lang_dat$juv == "adult-adult"])
t.test(lang_dat$fst[lang_dat$juv == "juv-juv"],lang_dat$fst[lang_dat$juv == "juv-adult"])
t.test(lang_dat$fst[lang_dat$juv == "juv-adult"],lang_dat$fst[lang_dat$juv == "adult-adult"])

lmar_dat <- lates_fst_ind_all %>%
  filter(!is.na(fst) & sampling_loc.pop1 == sampling_loc.pop2) %>%
  filter(spp == "lmar")

t.test(lmar_dat$fst[lmar_dat$juv == "juv-juv"],lmar_dat$fst[lmar_dat$juv == "adult-adult"])
t.test(lmar_dat$fst[lmar_dat$juv == "juv-juv"],lmar_dat$fst[lmar_dat$juv == "juv-adult"])
t.test(lmar_dat$fst[lmar_dat$juv == "juv-adult"],lmar_dat$fst[lmar_dat$juv == "adult-adult"])

lmic_dat <- lates_fst_ind_all %>%
  filter(!is.na(fst) & sampling_loc.pop1 == sampling_loc.pop2) %>%
  filter(spp == "lmic")

t.test(lmic_dat$fst[lmic_dat$juv == "juv-juv"],lmic_dat$fst[lmic_dat$juv == "adult-adult"])
t.test(lmic_dat$fst[lmic_dat$juv == "juv-juv"],lmic_dat$fst[lmic_dat$juv == "juv-adult"])
t.test(lmic_dat$fst[lmic_dat$juv == "juv-adult"],lmic_dat$fst[lmic_dat$juv == "adult-adult"])

lates_fst_ind_all %>% 
  filter(!is.na(fst) & sampling_loc.pop1 == sampling_loc.pop2) %>%
  ggplot() +
    geom_jitter(aes(x=juv,y=fst,fill=juv),size=4,pch=21, alpha=0.7,width=0.1,height=0) +
    geom_boxplot(aes(x=juv,y=fst),fill=NA,outlier.color=NA) +
    theme_custom() +
    scale_fill_npg() +
    facet_wrap(~spp,scales="free")


## Q5. Are juveniles collected at the same time of year more closely related?
lates_fst_ind_all %>% 
  filter(!is.na(fst) & juv == "juv-juv" & sampling_season.pop1 != "" &
           sampling_season.pop2 != "") %>%
  ggplot() +
  geom_jitter(aes(x=combined_season,y=fst, fill=juv),size=4,pch=21, alpha=0.7,width=0.1,height=0) +
  geom_boxplot(aes(x=combined_season,y=fst),fill=NA,outlier.color=NA) +
  theme_custom() +
  scale_fill_npg() +
  facet_wrap(~spp,scales="free")

lates_fst_ind_all %>% 
  filter(!is.na(fst) & juv == "juv-juv" & sampling_season.pop1 != "" &
           sampling_season.pop2 != "") %>%
  mutate(same_season = case_when(sampling_season.pop1 == sampling_season.pop2 ~ "same",
                                 TRUE ~ "different")) %>%
  ggplot() +
  geom_jitter(aes(x=same_season,y=fst, fill=juv),size=4,pch=21, alpha=0.7,width=0.1,height=0) +
  geom_boxplot(aes(x=same_season,y=fst),fill=NA,outlier.color=NA) +
  theme_custom() +
  scale_fill_npg() +
  facet_wrap(~spp,scales="free")

lmic_dat2 <- lates_fst_ind_all %>%
  filter(!is.na(fst) & juv == "juv-juv" & sampling_season.pop1 != "" &
           sampling_season.pop2 != "") %>%
  mutate(same_season = case_when(sampling_season.pop1 == sampling_season.pop2 ~ "same",
                                 TRUE ~ "different")) 

t.test(lmic_dat2$fst[lmic_dat2$same_season == "same"],lmic_dat2$fst[lmic_dat2$same_season == "different"])
