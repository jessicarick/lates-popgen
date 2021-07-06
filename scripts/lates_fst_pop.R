###################
## Script for fst calculations
##
## Analysis for Rick et al., in prep
## Written by J. Rick, jrick@uwyo.edu
## Last update: Summer 2021
###################

##################### 
## Loading packages and data
#####################
## Load required packages
source("packages_funcs.R")

## Importing the SNP data and metadata
## working with the dataset created with the following filters (from VCFtools): 
## missing data allowed = 50%, MAF = 0.01 cutoff

## Toggle for working with gbs vs rad data
type <- "gbs"
lates_vcfR <- read.vcfR("../data/lates_all_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lang_vcfR <-  read.vcfR("../data/lang_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmar_vcfR <-  read.vcfR("../data/lmar_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lmic_vcfR <-  read.vcfR("../data/lmic_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_vcfR <-  read.vcfR("../data/lsta_092320_0.5_maf0.01_thin90_dp5.recode.vcf")
lsta_rad_vcfR <- read.vcfR("../data/combined_lsta_noHets3-4_variants_0.5_maf0.01_dp5_thin90.recode.vcf")


#####################
## Cleaning up the data
#####################
vcfs <- c("lates_vcfR","lang_vcfR","lmar_vcfR","lmic_vcfR","lsta_vcfR","lsta_rad_vcfR")
for(v in c(lang_vcfR,lmar_vcfR,lmic_vcfR,lsta_vcfR)){
  latesgen <- vcfR2genlight(v)
  
  # clean up names and import associated metadata
    if (type == "gbs") {
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
  
  ## extracting genotype matrix from genlight
  # lates_alleles <- t(as.matrix(latesgen_nolowcov))
  # 
  # dim(lates_alleles)
  # head(colnames(lates_alleles)) ## to make sure that the names look good
  
  # #################
  # ## Initial PCA
  # #################
  # 
  # lates_pca <- do.pca(lates_alleles)
  # pcSummary <- summary(lates_pca)
  # scree <- plot(lates_pca, type="lines") 
  # 
  ## Pull metadata associated with individuals in the dataset
  pairedinfolates<-left_join(data.frame(Moran_FishID=indNames(latesgen_nolowcov)),
                             fishinfo,by="Moran_FishID",all.x=TRUE,all.y=FALSE)
  pairedinfolates[is.na(pairedinfolates)] <- "UNK"
  head(pairedinfolates) # make sure it paired okay


  # ## Combining fish info with PC results
  # if (type == "gbs") {
  #   pcaAll <- data.frame(names = pairedinfolates$Moran_FishID,
  #                        spp = factor(pairedinfolates$final_ID),
  #                        fieldID = factor(pairedinfolates$pheno_ID),
  #                        site = factor(pairedinfolates$sampling_loc),
  #                        library = factor(pairedinfolates$seq_library),
  #                        sex = factor(pairedinfolates$pheno_sex),
  #                        EV1 = lates_pca$x[,1],    # the first eigenvector
  #                        EV2 = lates_pca$x[,2],    # the second eigenvector
  #                        EV3 = lates_pca$x[,3],    # the third eigenvector
  #                        EV4 = lates_pca$x[,4],
  #                        EV5 = lates_pca$x[,5],
  #                        stringsAsFactors = FALSE)
  # } else if (type == "rad") {
  #   pcaAll <- data.frame(names = pairedinfolates$Indv_ID,
  #                        spp = factor(pairedinfolates$Spp),
  #                        site = factor(pairedinfolates$Location),
  #                        library = factor(pairedinfolates$Library),
  #                        sex = factor(pairedinfolates$Sex),
  #                        EV1 = lates_pca$x[,1],    # the first eigenvector
  #                        EV2 = lates_pca$x[,2],    # the second eigenvector
  #                        EV3 = lates_pca$x[,3],    # the third eigenvector
  #                        EV4 = lates_pca$x[,4],
  #                        EV5 = lates_pca$x[,5],
  #                        stringsAsFactors = FALSE)
  # }
  # 
  # ## Now, plotting the PCA results
  # 
  # # colored by sampling site
  # par(oma=c(1,1,1,2), xpd=TRUE, mar=c(5.1, 4.1, 4.1, 1),mfrow=c(1,3))
  # plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$site],0.6),col=colors2[pcaAll$site],
  #      xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
  #      ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))
  # legend("topright",legend=levels(pcaAll$site),col=scales::alpha(colors2,0.6),border=NULL,pch=19,bty="n", cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)
  # 
  # plot(pcaAll$EV2, pcaAll$EV3, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$site],0.6),col=colors2[pcaAll$site],
  #      xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
  #      ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))
  # 
  # plot(pcaAll$EV3, pcaAll$EV4, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$site],0.6),col=colors2[pcaAll$site],
  #      xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
  #      ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))
  # 
  # # colored by library (to make sure library effects are gone)
  # plot.new()
  # par(mar=c(6,6,1,4),mfrow=c(1,3),xpd=TRUE, oma=c(2,2,1,2))
  # plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
  #      col=colors[pcaAll$library],
  #      xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
  #      ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
  #      cex.lab=3,cex.axis=2)
  # legend("bottomleft",legend=levels(pcaAll$library),col=scales::alpha(colors,0.6),border=NULL,pch=19,bty="n", cex=2, pt.cex=2, pt.lwd=2, horiz=FALSE)
  # 
  # plot(pcaAll$EV2, pcaAll$EV3, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
  #      col=colors[pcaAll$library],
  #      xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
  #      ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
  #      cex.lab=3,cex.axis=2)
  # 
  # plot(pcaAll$EV3, pcaAll$EV4, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
  #      col=colors[pcaAll$library],
  #      xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
  #      ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
  #      cex.lab=3,cex.axis=2)
  
  #####################
  ## Calculate Reich-Patterson FST between populations
  #####################
  
  pop(latesgen_nolowcov) <- pairedinfolates$sampling_loc
  pop(latesgen.nolowcov) <- pairedinfolates$sampling_loc
  
  # lates_FST_spp_preLib <- reich.fst(latesgen.nolowcov,
  #                                   bootstrap=100,
  #                                   plot=TRUE,
  #                                   verbose=TRUE)
  
  # lates_FST_spp <- reich.fst(latesgen_nolowcov,
  #                            bootstrap=100,
  #                            plot=TRUE,
  #                            verbose=TRUE)
  
  # for comparison
  #lates_fst_spp_dartR <- gl.fst.pop(latesgen_nolowcov)
  
  # calculate general stats by species
  lates_heterozyg <- gl.report.heterozygosity(latesgen_nolowcov,method="pop")
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
p <- lates_FST_all3 %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),
          min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lsta",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lmic",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lmar",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lang",10)) %>%
  add_row(pop1=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          pop2=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                 "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                 "Kasanga"),
          fst_estimate=rep(NA,10),min_CI=rep(NA,10),max_CI=rep(NA,10),
          spp=rep("Lsta_RAD",10)) %>%
  mutate(pop1 = factor(pop1,
                       levels=rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                "Kasanga","Cameroon","Congo","Dar"))),
         pop2 = factor(pop2,
                       levels=rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                "Kasanga","Cameroon","Congo","Dar"))),
         spp = factor(spp,levels=c("Lsta","Lsta_RAD","Lmic","Lmar","Lang")))  %>%
  filter(as.integer(pop1) <= as.integer(pop2) &
         spp != "Lsta_RAD") %>%
  # transform(spp = factor(spp,levels=c("Lsta","Lsta_RAD","Lmic","Lmar","Lang"))) %>%
ggplot(aes(x=pop1,y=pop2,fill=fst_estimate)) +
  geom_tile(color="black",size=0.1) +
  geom_text(aes(label = round(fst_estimate, 3)),
            size=4,family="Open Sans Light") +
  geom_text(data=hets.all.spp[hets.all.spp$spp != "Lsta_RAD",],aes(x=pop,y=pop,label=round(Ho,3)),
            inherit.aes=FALSE,col="white",size=4,family="Open Sans Light") +
  facet_wrap(~factor(spp,levels=c("Lsta","Lmic","Lmar","Lang")),nrow=1,drop=FALSE,strip.position="top") +
  theme_custom() +
  theme(axis.text.x = element_text(angle=45,hjust=0),
        strip.text = element_blank(),
        panel.spacing = unit(5, "lines"),
        panel.border = element_rect(color="white"),
        axis.line.x.top = element_line(color="black"),
        axis.line.y.left = element_line(color="black")) +
  labs(x=NULL,y=NULL) +
  #scale_fill_viridis(name="FST",option="plasma",direction=-1,end=0.6)
  scale_fill_gradient(low="#D6D3CC",high="#0baf8f") +
  scale_x_discrete(limits=c(rev(c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                  "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                  "Kasanga"))),
                   position="top") +
  scale_y_discrete(limits=c(rev(c("","Kagunga","Kigoma","N_Mahale","S_Mahale",
                                  "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                  "Kasanga"))))
    #scale_fill_brewer(palette="RdYlGn")

p + geom_text(data=totals,aes(y=11,x=pop,label=n,fill=NULL),col="red",family="Open Sans",size=5)

#------------------------------#
lmar_fst_byInd$fsts %>%
  as_data_frame() %>%
  mutate(pop1 = colnames(.)) %>%
  pivot_longer(cols=!starts_with("pop"),names_to="pop2",values_to="fst") %>%
  left_join(fishinfo,by=c("pop1" = "ind")) %>%
  left_join(fishinfo,by=c("pop2" = "ind"),suffix=c(".pop1",".pop2")) %>%
  mutate(diff_sl = abs(as.numeric(SL_mm.pop1)-as.numeric(SL_mm.pop2)),
         juv = case_when(juvenile.pop1 == "Y" & juvenile.pop2 == "Y" ~ "juv-juv",
                         juvenile.pop1 == "Y" & juvenile.pop2 == "N" ~ "juv-adult",
                         juvenile.pop1 == "N" & juvenile.pop2 == "Y" ~ "juv-adult",
                         juvenile.pop1 == "N" & juvenile.pop2 == "N" ~ "adult-adult",
                         juvenile.pop1 == "" | juvenile.pop2 == "" ~ "")) %>%
  filter(!is.na(fst) & juv != "") %>%
  #ggplot(aes(x=fst,fill=juv)) +
  ggdensity(x="fst",alpha=0.5,fill="juv",add="mean",rug=TRUE,palette="aaas")
#-------------------------------#


# convert to genlight object and define ploidy
latesgen <- vcfR2genlight(lates_vcfR)
ploidy(latesgen) <- 2
latesgen <- gl.compliance.check(latesgen)

# clean up names and import associated metadata
if (type == "gbs") {
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
  
  #pairedinfolates$Library[pairedinfolates$Library == "lates02lates03"] <- "lates03"
  
} else if (type == "rad") {
  col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_all/bwa_all_rad/aln_"))
  col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
  indNames(latesgen)<-col.names.clean
  
  fishinfo <- read.csv('lates_rad_metadata.csv',
                       header=TRUE,stringsAsFactors = FALSE)
  
  pairedinfolates <- left_join(data.frame(Moran_ID=col.names.clean),
                               fishinfo, by="Moran_ID", all.x=TRUE, all.y=FALSE)
  
  pop(latesgen) <- pairedinfolates$Library
}

# now, filter out any individuals with > 50% missing data
latesgen.nolowcov <- gl.filter.callrate(latesgen,method="ind",
                                        threshold=0.5,mono.rm=TRUE,
                                        recalc=FALSE,plot=TRUE,v=2) #TODO -- remove?

#####################
## Library effects filter
## filtering SNPs found only in one group or the other
#####################
pop(latesgen.nolowcov) <- pairedinfolates$seq_library

if (type == "gbs"){
  lates01 <- gl.keep.pop(latesgen.nolowcov,
                         "lates01", 
                         mono.rm = FALSE)
  lates01.highcalls <- gl.filter.callrate(lates01, 
                                          method="loc",
                                          threshold=0.75,
                                          recalc=FALSE,
                                          mono.rm = FALSE)
  lates03 <- gl.keep.pop(latesgen.nolowcov,
                         "lates03", mono.rm=FALSE)
  lates03.highcalls <- gl.filter.callrate(lates03,
                                          method="loc",
                                          threshold=0.75, 
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
                                         threshold=0.75,
                                         recalc=FALSE,
                                         mono.rm = FALSE)
  
  GQI132 <- gl.keep.pop(latesgen,
                        "GQI132")
  GQI132.highcalls <- gl.filter.callrate(GQI132, 
                                         method="loc",
                                         threshold=0.75,
                                         recalc=FALSE,
                                         mono.rm = FALSE)
  
  GQI133 <- gl.keep.pop(latesgen.nolowcov,
                        "GQI133", mono.rm=FALSE)
  GQI133.highcalls <- gl.filter.callrate(GQI133, 
                                         method="loc",
                                         threshold=0.75,
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

## extracting genotype matrix from genlight
lates_alleles <- t(as.matrix(latesgen_nolowcov))

dim(lates_alleles)
head(colnames(lates_alleles)) ## to make sure that the names look good

#################
## Initial PCA
#################

lates_pca <- do.pca(lates_alleles)
pcSummary <- summary(lates_pca)
scree <- plot(lates_pca, type="lines") 

## Pull metadata associated with individuals in the dataset
pairedinfolates<-left_join(data.frame(Moran_FishID=indNames(latesgen_nolowcov)),
                           fishinfo,by="Moran_FishID",all.x=TRUE,all.y=FALSE) 
pairedinfolates[is.na(pairedinfolates)] <- "UNK"
head(pairedinfolates) # make sure it paired okay


## Combining fish info with PC results 
if (type == "gbs") {
  pcaAll <- data.frame(names = pairedinfolates$Moran_FishID,
                       spp = factor(pairedinfolates$final_ID),
                       fieldID = factor(pairedinfolates$pheno_ID),
                       site = factor(pairedinfolates$sampling_loc),
                       library = factor(pairedinfolates$seq_library),
                       sex = factor(pairedinfolates$pheno_sex),
                       EV1 = lates_pca$x[,1],    # the first eigenvector
                       EV2 = lates_pca$x[,2],    # the second eigenvector
                       EV3 = lates_pca$x[,3],    # the third eigenvector
                       EV4 = lates_pca$x[,4],
                       EV5 = lates_pca$x[,5],
                       stringsAsFactors = FALSE)
} else if (type == "rad") {
  pcaAll <- data.frame(names = pairedinfolates$Indv_ID,
                       spp = factor(pairedinfolates$Spp),
                       site = factor(pairedinfolates$Location),
                       library = factor(pairedinfolates$Library),
                       sex = factor(pairedinfolates$Sex),
                       EV1 = lates_pca$x[,1],    # the first eigenvector
                       EV2 = lates_pca$x[,2],    # the second eigenvector
                       EV3 = lates_pca$x[,3],    # the third eigenvector
                       EV4 = lates_pca$x[,4],
                       EV5 = lates_pca$x[,5],
                       stringsAsFactors = FALSE)
}

## Now, plotting the PCA results

# colored by species
par(oma=c(1,1,1,2), xpd=TRUE, mar=c(5.1, 4.1, 4.1, 1),mfrow=c(1,3))
plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))
legend("topright",legend=levels(pcaAll$spp),col=scales::alpha(colors2,0.6),border=NULL,pch=19,bty="n", cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(pcaAll$EV2, pcaAll$EV3, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))

plot(pcaAll$EV3, pcaAll$EV4, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))

# colored by library (to make sure library effects are gone)
plot.new()
par(mar=c(6,6,1,4),mfrow=c(1,3),xpd=TRUE, oma=c(2,2,1,2))
plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
     col=colors[pcaAll$library],
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)
legend("bottomleft",legend=levels(pcaAll$library),col=scales::alpha(colors,0.6),border=NULL,pch=19,bty="n", cex=2, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(pcaAll$EV2, pcaAll$EV3, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
     col=colors[pcaAll$library],
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)

plot(pcaAll$EV3, pcaAll$EV4, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5), 
     col=colors[pcaAll$library],
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)

## 3D and interactive plot for investigating the data
library(scatterplot3d)
scatterplot3d(x=pcaAll$EV1,y=pcaAll$EV2,z=pcaAll$EV3,color=scales::alpha(c("#2a9d8f","#e76f51","#70567d","#e9c46a")[pcaAll$spp],0.4),pch=19,cex.symbols=2)

library(plotly)
p <- plot_ly(pcaAll, x = ~EV1, y = ~EV2, z = ~EV3, color = ~spp, colors = colors.vir,
             text = ~paste('ID:', names, '<br>Library:', library)) %>%
  add_markers() %>%
  
  layout(scene = list(xaxis = list(title = 
                                     paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
                                           sep="")),
                      yaxis = list(title = 
                                     paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
                                           sep="")),
                      zaxis = list(title = 
                                     paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", 
                                           sep=""))))
p

#####################
## Calculate Reich-Patterson FST between species
#####################

pop(latesgen_nolowcov) <- pairedinfolates$final_ID
pop(latesgen.nolowcov) <- pairedinfolates$final_ID

lates_FST_spp_preLib <- reich.fst(latesgen.nolowcov,
                           bootstrap=100,
                           plot=TRUE,
                           verbose=TRUE)

lates_FST_spp <- reich.fst(latesgen_nolowcov,
                           bootstrap=100,
                           plot=TRUE,
                           verbose=TRUE)

# for comparison
lates_fst_spp_dartR <- gl.fst.pop(latesgen_nolowcov)

# calculate general stats by species
lates_heterozyg <- gl.report.heterozygosity(latesgen_nolowcov,method="pop")





