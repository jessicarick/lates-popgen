#!/usr/bin Rscript

## script for plotting EMU pcas
library(tidyverse)
source(theme_custom.R)
library(RColorBrewer)
library(ggpubr)

spp_cols <- c("#2a9d8f","#e76f51","#70567d","#e9c46a")
loc_cols <- c("#854B47","#f25f5c","#f9a061",
              "#ffe066","#92ae83","#b0e8b3","#AAD3DA",
              "#2A8CB7","#833E88")


#args <- commandArgs(trailingOnly=TRUE)
dir <- "../../../results/emu_062021/"
args <- paste0(dir,c("lsta_092320_0.5_maf0.01_thin90_dp5",
          "lmic_092320_0.5_maf0.01_thin90_dp5",
          "lmar_092320_0.5_maf0.01_thin90_dp5",
          "lang_092320_0.5_maf0.01_thin90_dp5"))

## load information about missingness
## from vcftools --missing-indv option
lang.imiss <- read_table2("../data/lang_092320_0.5_maf0.01_thin90_dp5.imiss") %>%
  mutate(SPP = "Lang")
lmar.imiss <- read_table2("../data/lmar_092320_0.5_maf0.01_thin90_dp5.imiss") %>%
  mutate(SPP = "Lmar")
lmic.imiss <- read_table2("../data/lmic_092320_0.5_maf0.01_thin90_dp5.imiss") %>%
  mutate(SPP = "Lmic")
lsta.imiss <- read_table2("../data/lsta_092320_0.5_maf0.01_thin90_dp5.imiss") %>%
  mutate(SPP = "Lsta")
all.imiss <- rbind(lang.imiss,lmar.imiss,lmic.imiss,lsta.imiss)  %>%
  filter(F_MISS < 0.5)
ggdensity(all.imiss,x="F_MISS",fill="SPP",col="SPP",palette=spp_cols,rug=TRUE,add="mean") 

## load information about depth
## from vcftools --site-mean-depth option
lang.ldepth <- read_table2("../data/lang_092320_0.5_maf0.01_thin90_dp5.ldepth.mean") %>%
  mutate(SPP = "Lang")
lmar.ldepth <- read_table2("../data/lmar_092320_0.5_maf0.01_thin90_dp5.ldepth.mean") %>%
  mutate(SPP = "Lmar")
lmic.ldepth <- read_table2("../data/lmic_092320_0.5_maf0.01_thin90_dp5.ldepth.mean") %>%
  mutate(SPP = "Lmic")
lsta.ldepth <- read_table2("../data/lsta_092320_0.5_maf0.01_thin90_dp5.ldepth.mean") %>%
  mutate(SPP = "Lsta")
all.ldepth <- rbind(lang.ldepth,lmar.ldepth,lmic.ldepth,lsta.ldepth)
depth.stats <- all.ldepth %>%
  group_by(SPP) %>%
  summarise(depth.mean = mean(MEAN_DEPTH),
            min.depth = min(MEAN_DEPTH),
            ci.low = quantile(MEAN_DEPTH,0.025),
            ci.upp = quantile(MEAN_DEPTH,0.975))

ggdensity(all.ldepth,x="MEAN_DEPTH",fill="SPP",col="SPP",palette=spp_cols,rug=TRUE) +
  geom_segment(data=depth.stats,aes(x=depth.mean,xend=depth.mean,y=0,yend=0.075,col=SPP),lty=2,size=1.5) +
  geom_text(data=depth.stats,aes(x=depth.mean,y=0.08,col=SPP,label=round(depth.mean,1)),angle=90)

if (length(args) == 1) {
  base <- args[1]
  info <- read_csv("../data/lates_all_metadata.csv") %>%
    mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                                                   "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                                                   "Kasanga"))) %>%
    left_join(all.imiss,by=c("ind" = "FID")) 
  
  vec <- read_table2(paste0(base,".emu.eigenvecs"),col_names=FALSE) # Reads in eigenvectors
  val <- read_table2(paste0(base,".emu.eigenvals"),col_names="eigenval") %>%
    mutate(pct = eigenval/sum(eigenval))
  fam <- read_table2(paste0(base,".fam"),col_names=FALSE)
  colnames(fam)[1] <- "ind"
  
  combined <- fam %>%
    dplyr::select(ind) %>%
    left_join(info,by="ind") %>%
    bind_cols(vec)
  
  plot <- combined %>%
    ggplot(aes(x=PC1,y=PC2,col=sampling_loc)) +
    geom_point(size=8,alpha=0.5) + 
    theme_custom() +
    xlab(paste0("PC1 (",round(val$pct[1]*100,1),"%)")) +
    ylab(paste0("PC2 (",round(val$pct[2]*100,1),"%)")) +
    scale_colour_manual(values=colors,na.translate=FALSE)
  
  ggsave(paste(base,"_emu_pca.png"),plot=plot)
  
} else if (length(args) > 1) {
  plotlist <-  vector('list', length(args))
  miss.plotlist <- vector('list',2*length(args))
  year.plotlist <- vector('list',length(args))
  
  info <- read_csv("../data/lates_all_metadata.csv") %>%
    mutate(sampling_loc = case_when(sampling_loc == "N_Mahale" ~ "North Mahale",
                                    sampling_loc == "S_Mahale" ~ "South Mahale",
                                    TRUE ~ sampling_loc)) %>%
    mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","North Mahale","South Mahale",
                                                       "Ikola","Mpinbwe","Kirando","Wampembe",
                                                       "Kasanga"))) %>%
    left_join(all.imiss,by=c("ind" = "FID"))
  
  for (i in 1:length(args)){
    base <- args[i]
    spp <- case_when(grepl("lang",base) ~ "L. angustifrons",
                     grepl("lmar",base) ~ "L. mariae",
                     grepl("lmic",base) ~ "L. microlepis",
                     grepl("lsta",base) ~ "L. stappersii")
    
    
    vec <- read_table2(paste0(base,".emu.eigenvecs"),col_names=FALSE) # Reads in eigenvectors
    val <- read_table2(paste0(base,".emu.eigenvals"),col_names="eigenval") %>%
      mutate(pct = eigenval/sum(eigenval))
    fam <- read_table2(paste0(base,".fam"),col_names=FALSE)
    colnames(fam)[1] <- "ind"
    
    combined <- fam %>%
      dplyr::select(ind) %>%
      left_join(info,by="ind") %>%
      bind_cols(vec) %>%
      drop_na(all_of(c("sampling_loc","final_ID")))
    
    plot <- combined %>%
      ggplot(aes(x=X1,y=X2,col=sampling_loc,fill=sampling_loc)) +
      geom_point(pch=21,size=7,alpha=0.8) + 
      theme_custom() +
      xlab(paste0("PC1 (",round(val$pct[1]*100,2),"%)")) +
      ylab(paste0("PC2 (",round(val$pct[2]*100,2),"%)")) +
      #scale_color_discrete(drop=FALSE) 
      # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
      #          size=10)
      scale_colour_discrete(drop=FALSE) +
      scale_fill_discrete(drop=FALSE) + 
      guides(colour = guide_legend(nrow = 2)) +
      theme(legend.text = element_text(size=18),
            legend.position = "bottom",
            legend.margin = margin(15, 0, 15, 0),
            plot.margin = margin(30,5,10,5))
    
    plot.year <- combined %>%
      ggplot(aes(x=X1,y=X2,col=factor(sampling_year,
                                         levels=c("2019", "2018", "2017", "2016", "2015", "2012", "2011", "2010", "2002", "2001")))) +
      geom_point(size=7,alpha=0.8,stroke=2,
                 aes(shape=factor(sampling_year,
                                  levels=c("2019", "2018", "2017", "2016", "2015", "2012", "2011", "2010", "2002", "2001")))) + 
      theme_custom() +
      xlab(paste0("PC1 (",round(val$pct[1]*100,2),"%)")) +
      ylab(paste0("PC2 (",round(val$pct[2]*100,2),"%)")) +
      #scale_color_discrete(drop=FALSE) 
      # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
      #          size=10)
      scale_colour_discrete(drop=FALSE) +
      scale_fill_discrete(drop=FALSE) + 
      #guides(colour = guide_legend(nrow = 2)) +
      theme(legend.text = element_text(size=18),
            legend.position = "right",
            legend.margin = margin(15, 0, 15, 0),
            plot.margin = margin(30,5,10,5)) +
      scale_shape_manual(values=c(1,2,3,4,5,7,15,16,17,18),
                         drop=FALSE)
    
    plot.miss <- combined %>%
      ggplot(aes(x=X1,y=F_MISS,col=sampling_loc,fill=sampling_loc)) +
      geom_point(pch=21,size=7,alpha=0.8) + 
      theme_custom() +
      xlab(paste0("PC1 (",round(val$pct[1]*100,2),"%)")) +
      ylab(paste0("Missingness (%)")) +
      #scale_color_discrete(drop=FALSE) 
      # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
      #          size=10)
      scale_colour_discrete(drop=FALSE) +
      scale_fill_discrete(drop=FALSE) + 
      guides(colour = guide_legend(nrow = 2)) +
      theme(legend.text = element_text(size=18),
            legend.position = "bottom",
            legend.margin = margin(15, 0, 15, 0),
            plot.margin = margin(30,5,10,5))
    plot.miss2 <- combined %>%
      ggplot(aes(x=X2,y=F_MISS,col=sampling_loc,fill=sampling_loc)) +
      geom_point(pch=21,size=7,alpha=0.8) + 
      theme_custom() +
      xlab(paste0("PC2 (",round(val$pct[2]*100,2),"%)")) +
      ylab(paste0("Missingness (%)")) +
      #scale_color_discrete(drop=FALSE) 
      # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
      #          size=10)
      scale_colour_discrete(drop=FALSE) +
      scale_fill_discrete(drop=FALSE) + 
      guides(colour = guide_legend(nrow = 2)) +
      theme(legend.text = element_text(size=18),
            legend.position = "bottom",
            legend.margin = margin(15, 0, 15, 0),
            plot.margin = margin(30,5,10,5))
    
    plotlist[[i]] <- plot
    miss.plotlist[[2*i-1]] <- plot.miss
    miss.plotlist[[2*i]] <- plot.miss2
    year.plotlist[[i]] <- plot.year
  }
  
  
  plot_together <- ggarrange(plotlist=plotlist,nrow=1,
                             legend="bottom",common.legend = TRUE,
                             labels=NULL)
  plot_together_miss <- ggarrange(plotlist=miss.plotlist,nrow=2,ncol=4,
                                  legend="bottom",common.legend = TRUE,
                                  labels = c("(A) L. stappersii","","(B) L. microlepis","",
                                             "(C) L. mariae","","(D) L. angustifrons"),
                                  font.label=list(family="Open Sans Light",size=18))
  plot_together_year <- ggarrange(plotlist=year.plotlist,legend="right",
                                  common.legend=TRUE,
                                  labels = c("(A) L. stappersii","(B) L. microlepis",
                                             "(C) L. mariae","(D) L. angustifrons"),
                                  font.label=list(family="Open Sans Light",size=18))
  #ggsave("emu_pca.png",plot=plot_together)
  
  png("emu_pca.png",width=1500,height=375)
  print(plot_together)
  dev.off()  
}


#########################
## make a plot for the species combined dataset
#########################

base <- "../../../lates_popgen/results/emu_062021/lates_all_092320_0.5_maf0.01_thin90_dp5"
info <- read_csv("../../../lates_popgen/scripts/lates_popgen_github/data/lates_all_metadata.csv") %>%
  mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                                     "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                                     "Kasanga")))

vec <- read_table2(paste0(base,".emu.eigenvecs"),col_names=paste0("PC",seq(1,150))) # Reads in eigenvectors
val <- read_table2(paste0(base,".emu.eigenvals"),col_names="eigenval") %>%
  mutate(pct = eigenval/sum(eigenval))
fam <- read_table2(paste0(base,".fam"),col_names=FALSE)
colnames(fam)[1] <- "ind"

combined <- fam %>%
  dplyr::select(ind) %>%
  left_join(info,by="ind") %>%
  bind_cols(vec)

plot <- combined %>%
  ggplot(aes(x=PC1,y=PC2,col=final_ID)) +
  geom_point(size=8,alpha=0.5) + 
  theme_custom() +
  xlab(paste0("PC1 (",round(val$pct[1]*100,1),"%)")) +
  ylab(paste0("PC2 (",round(val$pct[2]*100,1),"%)")) +
  scale_colour_manual(values=spp_cols,na.translate=FALSE)
plot2 <- combined %>%
  ggplot(aes(x=PC2,y=PC3,col=final_ID)) +
  geom_point(size=8,alpha=0.5) + 
  theme_custom() +
  xlab(paste0("PC2 (",round(val$pct[2]*100,1),"%)")) +
  ylab(paste0("PC3 (",round(val$pct[3]*100,1),"%)")) +
  scale_colour_manual(values=spp_cols,na.translate=FALSE)
plot3 <- combined %>%
  ggplot(aes(x=PC3,y=PC4,col=final_ID)) +
  geom_point(size=8,alpha=0.5) + 
  theme_custom() +
  xlab(paste0("PC3 (",round(val$pct[3]*100,1),"%)")) +
  ylab(paste0("PC4 (",round(val$pct[4]*100,1),"%)")) +
  scale_colour_manual(values=spp_cols,na.translate=FALSE)

plots.combined <- ggarrange(plot,plot2,plot3,nrow=1,common.legend = TRUE)
ggsave(paste(base,"_emu_pca.png"),plot=plot)
ggsave(paste(base,"_plots.combined.png"),plot=plots.combined,
       width=15,height=5,units="in")

       
#-------------#
ggarrange(plotlist[[1]],lsta3,
          plotlist[[2]],lmic3,
          plotlist[[3]],lmar3,
          plotlist[[4]],lang3,
          common.legend=TRUE,
          ncol=2,nrow=4,align="v",
          legend="bottom") + guides(colour = guide_legend(nrow = 2))

#------------#

## lates rad emu plot
vec <- read_table2(paste0(dir,"combined_lsta_noHets3-4_variants_0.5_maf0.01.emu.eigenvecs"),col_names=FALSE) # Reads in eigenvectors
val <- read_table2(paste0(dir,"combined_lsta_noHets3-4_variants_0.5_maf0.01.emu.eigenvals"),col_names="eigenval") %>%
  mutate(pct = eigenval/sum(eigenval))
fam <- read_table2(paste0(dir,"combined_lsta_noHets3-4_variants_0.5_maf0.01_dp5_thin90.fam"),col_names=FALSE)
colnames(fam)[1] <- "ind"

info <- read_csv("../data/lates_all_metadata.csv") %>%
  mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                                     "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                                     "Kasanga")))

combined <- fam %>%
  dplyr::select(ind) %>%
  left_join(info,by="ind") %>%
  bind_cols(vec) %>%
  drop_na(all_of(c("sampling_loc","final_ID"))) %>%
  filter(ind != "KAG22.GQI132" & ind != "KAG22.GQI133")

plot1 <- combined %>%
  ggplot(aes(x=X1,y=X2,col=sampling_loc,fill=sampling_loc)) +
  geom_point(pch=21,size=7,alpha=0.8) + 
  theme_custom() +
  xlab(paste0("PC1 (",round(val$pct[1]*100,2),"%)")) +
  ylab(paste0("PC2 (",round(val$pct[2]*100,2),"%)")) +
  #scale_color_discrete(drop=FALSE) 
  # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
  #          size=10)
  scale_colour_discrete(drop=FALSE) +
  scale_fill_discrete(drop=FALSE) + 
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.text = element_text(size=18),
        legend.position = "bottom",
        legend.margin = margin(15, 0, 15, 0),
        plot.margin = margin(30,5,10,5))

plot2 <- combined %>%
  ggplot(aes(x=X2,y=X3,col=sampling_loc,fill=sampling_loc)) +
  geom_point(pch=21,size=7,alpha=0.8) + 
  theme_custom() +
  xlab(paste0("PC2 (",round(val$pct[2]*100,2),"%)")) +
  ylab(paste0("PC3 (",round(val$pct[3]*100,2),"%)")) +
  #scale_color_discrete(drop=FALSE) 
  # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
  #          size=10)
  scale_colour_discrete(drop=FALSE) +
  scale_fill_discrete(drop=FALSE) + 
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.text = element_text(size=18),
        legend.position = "bottom",
        legend.margin = margin(15, 0, 15, 0),
        plot.margin = margin(30,5,10,5))

plot3 <- combined %>%
  ggplot(aes(x=X3,y=X4,col=sampling_loc,fill=sampling_loc)) +
  geom_point(pch=21,size=7,alpha=0.8) + 
  theme_custom() +
  xlab(paste0("PC3 (",round(val$pct[3]*100,2),"%)")) +
  ylab(paste0("PC4 (",round(val$pct[4]*100,2),"%)")) +
  #scale_color_discrete(drop=FALSE) 
  # annotate("text",y=max(combined$PC2),x=max(combined$PC1),label=spp,hjust=1,
  #          size=10)
  scale_colour_discrete(drop=FALSE) +
  scale_fill_discrete(drop=FALSE) + 
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.text = element_text(size=18),
        legend.position = "bottom",
        legend.margin = margin(15, 0, 15, 0),
        plot.margin = margin(30,5,10,5))

ggarrange(plot1,plot2,plot3,nrow=1,common.legend=TRUE)
