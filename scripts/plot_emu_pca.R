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

args <- paste0(dir,c("lsta_092320_0.5_maf0.01_thin90_dp5",
          "lmic_092320_0.5_maf0.01_thin90_dp5",
          "lmar_092320_0.5_maf0.01_thin90_dp5",
          "lang_092320_0.5_maf0.01_thin90_dp5"))

if (length(args) == 1) {
  base <- args[1]
  info <- read_csv("../../../lates_popgen/scripts/lates_popgen_github/data/lates_all_metadata.csv") %>%
    mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                                                   "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                                                   "Kasanga","Cameroon","Congo","Dar")))
  
  vec <- read_table2(paste0(base,".emu.eigenvecs"),col_names=paste0("PC",seq(1,20))) # Reads in eigenvectors
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
  
  info <- read_csv("../../../lates_popgen/scripts/lates_popgen_github/data/lates_all_metadata.csv") %>%
    mutate(sampling_loc = case_when(sampling_loc == "N_Mahale" ~ "North Mahale",
                                    sampling_loc == "S_Mahale" ~ "South Mahale",
                                    TRUE ~ sampling_loc)) %>%
    mutate(sampling_loc = factor(sampling_loc,levels=c("Kagunga","Kigoma","North Mahale","South Mahale",
                                                       "Ikola","Mpinbwe","Kirando","Wampembe",
                                                       "Kasanga")))
  
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
      guides(colour = guide_legend(nrow = 1)) +
      theme(legend.text = element_text(size=18),
            legend.position = "bottom",
            legend.margin=margin(15, 0, 15, 0))
    
    plotlist[[i]] <- plot
  }
  
  
  plot_together <- ggarrange(plotlist=plotlist,nrow=1,
                             legend="bottom",common.legend = TRUE,
                             labels=NULL)
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
                                                     "Kasanga","Cameroon","Congo","Dar")))

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

       