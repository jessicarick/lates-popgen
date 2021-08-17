###################
## Script for fst calculations
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

## Importing the SNP data and metadata
## working with the datasets created with the following filters (from VCFtools): 
## missing data allowed = 50%, MAF = 0.01 cutoff, min depth = 5, thinned 90bp

lates_vcfR <- read.vcfR("../../lates_popgen_github_copy/data/lates_all_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lang_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lang_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmar_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lmar_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmic_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lmic_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_vcfR <-  read.vcfR("../../lates_popgen_github_copy/data/lsta_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_rad_vcfR <- read.vcfR("../data/combined_lsta_noHets3-4_variants_0.5_maf0.01_dp5_thin90.recode.vcf")

#####################
## Cleaning up the data
#####################
vcfs <- c(lates_vcfR,lang_vcfR,lmar_vcfR,lmic_vcfR,lsta_vcfR,lsta_rad_vcfR)
vcf_names <- c("lates_vcfR","lang_vcfR","lmar_vcfR","lmic_vcfR","lsta_vcfR","lsta_rad_vcfR")
for(i in 1:length(vcfs)){
  v <- vcfs[i]
  latesgen <- vcfR2genlight(v)
  type <- case_when(vcf_names[i] == "lsta_rad_vcfR" ~ "rad",
                    TRUE ~ "gbs")
  
    # clean up names and import associated metadata
    if (type == "gbs") {
    col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
    col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
    col.names.clean[col.names.clean == "CEW16_135_2"] <- "CEW16_135"
    indNames(latesgen)<-col.names.clean
    
    fishinfo <- read.csv('../data/lates_all_metadata.csv',
                         header=TRUE,stringsAsFactors = FALSE) %>%
      mutate(Moran_FishID = ind,
             SL_mm = as.numeric(SL_mm),
             TL_mm = as.numeric(TL_mm))
    
    pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                                 fishinfo[,-1], by="Moran_FishID",
                                 all.x=TRUE, all.y=FALSE)
    
    #pairedinfolates$Library[pairedinfolates$Library == "lates02lates03"] <- "lates03"
    pop(latesgen) <- pairedinfolates$seq_library
    
  } else if (type == "rad") {
    col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_all/bwa_all_rad/aln_"))
    col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
    indNames(latesgen)<-col.names.clean
    
    fishinfo <- read.csv('lates_rad_metadata.csv',
                         header=TRUE,stringsAsFactors = FALSE)
    
    pairedinfolates <- left_join(data.frame(Moran_ID=col.names.clean),
                                 fishinfo, by="Moran_ID", all.x=TRUE, all.y=FALSE)
    
    pop(latesgen) <- pairedinfolates$seq_library
  }
  
  
  # convert to genlight object and define ploidy
  # then filter out any individuals with > 50% missing data
  ploidy(latesgen) <- 2
  latesgen <- gl.compliance.check(latesgen)
  latesgen.nolowcov <- gl.filter.callrate(latesgen,method="ind",
                                          threshold=0.5,mono.rm=TRUE,
                                          recalc=FALSE,plot=TRUE,v=2) #TODO -- remove?
  
  #####################
  ## Library effects filter
  ## filtering SNPs found only in one group or the other
  #####################
  
  if (type == "gbs"){
    lates01 <- gl.keep.pop(latesgen.nolowcov,
                           "lates01", 
                           mono.rm = FALSE)
    lates01.highcalls <- gl.filter.callrate(lates01, 
                                            method="loc",
                                            threshold=0.5,
                                            recalc=FALSE,
                                            mono.rm = FALSE)
    lates03 <- gl.keep.pop(latesgen.nolowcov,
                           "lates03", mono.rm=FALSE)
    lates03.highcalls <- gl.filter.callrate(lates03,
                                            method="loc",
                                            threshold=0.5, 
                                            recalc=FALSE, 
                                            mono.rm=FALSE)
    
    snps.all <- Reduce(intersect, 
                       list(locNames(lates01.highcalls),
                            locNames(lates03.highcalls)))
    
  } else if (type == "rad") {
    
    GQI128 <- gl.keep.pop(latesgen.nolowcov,
                          "GQI128", mono.rm = FALSE)
    GQI128.highcalls <- gl.filter.callrate(GQI128, 
                                           method="loc",
                                           threshold=0.5,
                                           recalc=FALSE,
                                           mono.rm = FALSE)
    
    GQI132 <- gl.keep.pop(latesgen,
                          "GQI132")
    GQI132.highcalls <- gl.filter.callrate(GQI132, 
                                           method="loc",
                                           threshold=0.5,
                                           recalc=FALSE,
                                           mono.rm = FALSE)
    
    GQI133 <- gl.keep.pop(latesgen.nolowcov,
                          "GQI133", mono.rm=FALSE)
    GQI133.highcalls <- gl.filter.callrate(GQI133, 
                                           method="loc",
                                           threshold=0.5,
                                           recalc=FALSE,
                                           mono.rm = FALSE)
    
    snps.all <- Reduce(intersect, list(locNames(GQI128.highcalls),
                                       locNames(GQI132.highcalls),
                                       locNames(GQI133.highcalls)))
  }
  
  latesgen.all <- gl.keep.loc(latesgen.nolowcov, snps.all)
  
  # remove any indv that now have >50% missing data
  latesgen_nolowcov <- gl.filter.callrate(latesgen.all,
                                          method="ind",
                                          threshold=0.5,
                                          mono.rm=T,
                                          recalc=FALSE,
                                          plot=TRUE,v=2)
  

  ## Pull metadata associated with individuals in the dataset
  pairedinfolates <- left_join(data.frame(Moran_FishID=indNames(latesgen_nolowcov)),
                             fishinfo,by="Moran_FishID",all.x=TRUE,all.y=FALSE)
  pairedinfolates[is.na(pairedinfolates)] <- "UNK"
  head(pairedinfolates) # make sure it paired okay

  
  #####################
  ## Calculate site-level heterozygosity
  #####################

  pop(latesgen_nolowcov) <- pairedinfolates$sampling_loc
  spp <- pairedinfolates$final_ID[1]
  lates_heterozyg <- gl.report.heterozygosity(latesgen_nolowcov,method="pop") %>%
    mutate(spp=spp)
  assign(paste0(spp,"_heterozyg"),lates_heterozyg)
  plot(Ho~He,data=lates_heterozyg)
  abline(a=0,b=1,lty=2)
}

lates_FST_spp$bootstraps$spp <- "Lsta_RAD"
lates_FST_all_tmp <- rbind(lates_FST_all,lates_FST_spp$bootstraps)
lates_FST_all <- lates_FST_all_tmp

lates_FST_all$fst_estimate <- as.numeric(lates_FST_all$fst_estimate)

lates_FST_all2 <- lates_FST_all[,c(2,1,3:5,106)]
colnames(lates_FST_all2) <- c("pop1","pop2","fst_estimate","min_CI","max_CI","spp")

lates_FST_all3 <- rbind(lates_FST_all[,c(1:5,106)],lates_FST_all2)
lates_FST_all3$fst_estimate[lates_FST_all3$fst_estimate < 0] <- 0

################################
## write file with fst estimates
write.csv(lates_FST_all3,file="../data/lates_fst_by_sampling_site.csv",row.names=FALSE)
################################

totals <- fishinfo %>%
  group_by(final_ID,sampling_loc) %>%
  tally() %>%
  filter(final_ID != "#N/A") %>%
  mutate(spp = factor(final_ID, levels=c("Lsta","Lmar","Lmic","Lang")),
         pop = factor(sampling_loc, levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                             "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                             "Kasanga","Total")))

lates_FST_all4 <- lates_FST_all3 %>%
  left_join(totals,by=c("pop1" = "sampling_loc", "spp" = "final_ID")) %>%
  left_join(totals,by=c("pop2" = "sampling_loc", "spp" = "final_ID"), suffix=c("pop1","pop2")) %>%
  mutate(ntotal = npop1 + npop2)

################################
## plot heatmap of fst values by species
################################
hets.all.spp <- Lsta_heterozyg %>%
  bind_rows(Lang_heterozyg) %>%
  bind_rows(Lmar_heterozyg) %>%
  bind_rows(Lmic_heterozyg)

lates_FST_all3 <- read_csv("../data/lates_fst_by_sampling_site.csv")
lates_FST_plot_data <- lates_FST_all3 %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),
          min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lsta",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lmic",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lmar",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lang",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lsta_RAD",10)) %>%
  mutate(pop1 = factor(pop1,
                       levels=rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                "Ikola","Mpinbwe","Kirando","Wampembe",
                                "Kasanga","Cameroon","Congo","Dar"))),
         pop2 = factor(pop2,
                       levels=rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                "Ikola","Mpinbwe","Kirando","Wampembe",
                                "Kasanga","Cameroon","Congo","Dar"))),
         spp = factor(spp,levels=c("Lsta","Lsta_RAD","Lmic","Lmar","Lang")))  %>%
  filter(as.integer(pop1) <= as.integer(pop2) &
         spp != "Lsta_RAD") %>%
  # transform(spp = factor(spp,levels=c("Lsta","Lsta_RAD","Lmic","Lmar","Lang"))) %>%
  mutate(fill_col=case_when(min_CI <= 0 ~ NA_real_,
                            fst_estimate > 0 ~ fst_estimate))

p <- ggplot(data=lates_FST_plot_data) +
  geom_tile(aes(x=pop1,y=pop2,fill=fill_col),color="black",size=0.1) +
  geom_text(aes(x=pop1,y=pop2,label = round(fst_estimate, 3)),
            size=4,family="Open Sans") +
  geom_tile(data=filter(lates_FST_plot_data,pop1 == pop2),aes(x=pop1,y=pop2),color="black",fill="gray30",size=0.1) +
  geom_text(data=hets.all.spp[hets.all.spp$spp != "Lsta_RAD",],aes(x=pop,y=pop,label=round(Ho,3)),
            inherit.aes=FALSE,col="white",size=4,family="Open Sans") +
  facet_wrap(~factor(spp,levels=c("Lsta","Lmic","Lmar","Lang")),nrow=2,drop=FALSE,strip.position="top",scales="free") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45,hjust=0),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.border = element_rect(color="white"),
        axis.line.x.top = element_line(color="black"),
        axis.line.y.left = element_line(color="black"),
        legend.position="none",
        plot.margin = unit(c(0,1,0,0), "cm")) +
  labs(x=NULL,y=NULL) +
  # scale_fill_viridis(name="FST",option="inferno",begin=0.6,direction=-1,na.value="#D6D3CC",
  #                    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_fill_gradient(low="#f6e551",high="#e76f51",na.value="#D6D3CC",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_discrete(limits=c(rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                  "Ikola","Mpinbwe","Kirando","Wampembe",
                                  "Kasanga"))),
                   position="top") +
  scale_y_discrete(limits=c(rev(c("","Kagunga","Kigoma","N_Mahale","S_Mahale",
                                  "Ikola","Mpinbwe","Kirando","Wampembe",
                                  "Kasanga"))))
    #scale_fill_brewer(palette="RdYlGn")

p + geom_text(data=totals,aes(y=10,x=pop,label=n,fill=NULL),col="red",family="Open Sans",size=5)


###################################
#---------------------------------#
## Calculate FST between species ##
#---------------------------------#
###################################
  
# convert to genlight object and define ploidy
latesgen <- vcfR2genlight(lates_vcfR)
ploidy(latesgen) <- 2
latesgen <- gl.compliance.check(latesgen)

# clean up names and import associated metadata
col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
col.names.clean[col.names.clean == "CEW16_135_2"] <- "CEW16_135"
indNames(latesgen)<-col.names.clean

fishinfo <- read.csv('../data/lates_all_metadata.csv',
                     header=TRUE,stringsAsFactors = FALSE) %>%
  mutate(Moran_FishID = ind)
  
pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                             fishinfo[,-1], by="Moran_FishID",
                             all.x=TRUE, all.y=FALSE)

## Now, calculate Reich-Patterson FST between species
pop(latesgen) <- pairedinfolates$final_ID

lates_FST_spp <- reich.fst(latesgen,
                           bootstrap=100,
                           plot=TRUE,
                           verbose=TRUE)

# for comparison
lates_fst_spp_dartR <- gl.fst.pop(latesgen_nolowcov, nboots=100)

# calculate general stats by species
lates_heterozyg <- gl.report.heterozygosity(latesgen_nolowcov,method="pop")





