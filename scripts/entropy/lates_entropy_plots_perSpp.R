#####################
## Script for plotting entropy results
## for individual species, k=2 to k=8
##
## Written by J.rick, jrick@uwyo.edu
## Last updated Fall 2020
#####################

library(rhdf5)
library(abind)
library(RColorBrewer)
library(tidyverse)
library(memisc)
colors <- brewer.pal(3, "Set3")

## Accompanying metadata, change for your project
details <- read_csv("../data/lates_all_checkHets_admix_comparison.csv",
                    na=c("","#N/A"))[,1:16]
colnames(details)[1] <- "names"

for (spp in c("lsta","lmic","lmar","lang")){

for (i in 2:8){
  k <- i
  
  names <- read.table(paste0("../data/admixture/",spp,".ind"), header=F)
  colnames(names) <- c("names")
  fishinfo <- left_join(names, details, by.x="names", by.y="names",
                        all.x=T, all.y=F)
  
  ## Extract relevant parameters from .hdf5 format
  data1.q <- h5read(paste("../data/admixture/",spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep1.hdf5",sep=""), "q")
  data2.q <- h5read(paste("../data/admixture/",spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep1.hdf5",sep=""), "q")
  data3.q <- h5read(paste("../data/admixture/",spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep1.hdf5",sep=""), "q")
  
  ## Combine chains for each parameter estimated, take mean
  allq <- abind(data1.q,data2.q,data3.q,along=1)
  q <- t(apply(allq, 2:3, mean))
  
  ## Get credible intervals and calculate mean
  q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
  q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
  print(mean(q.ci.width))
  print(median(q.ci.width))
  
  ## Plot point estimates for q 
  # pdf(paste0("q_hist_",spp,"k",k,".pdf"), width=6.5, height=9)
  # par(mfrow=c(2,2))
  # for(i in 1:length(unique(fishinfo$Field.ID))){
  #   hist(q[which(fishinfo$Field.ID ==
  #                unique(fishinfo$Field.ID)[i]),1], 
  #        xlab="proportion of ancestry (q)", 
  #        main=unique(fishinfo$Field.ID)[i],
  #        col="gray90", xlim=c(0,1))
  # }
  # dev.off()
  
  # pdf(paste0(spp,"_conv_plots_k",k,".pdf"))
  # par(mfrow=c(4,4))
  # for (i in 1:nrow(names)){
  #   plot(data1.q[,1,i],type="l",main=names[i,1],ylim=c(0,1))
  #   lines(data2.q[,1,i],col="blue")
  #   lines(data3.q[,1,i],col="red")
  # }
  # dev.off()
  
  
  # stacked barplot
  q.df <- data.frame(names=fishinfo$names,
                     site=factor(as.character(fishinfo$sampling_loc), 
                                 levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar")),
                     spp=factor(fishinfo$final_ID),
                     #lib=factor(fishinfo$Library),
                     #plate=factor(fishinfo$Plate),
                     #dup=fishinfo$duplicate,
                     #sex=fishinfo$Sex,
                     q[,1:k]
  )
  
  # q.df$assign <- factor(case_when(q.df$X1 > 0.6 ~ "q1",
  #                          q.df$X2 > 0.6 ~ "q2",
  #                          q.df$X3 > 0.6 ~ "q3",
  #                          q.df$X4 > 0.6 ~ "q4",
  #                          q.df$X5 > 0.6 ~ "q5",
  #                          TRUE ~ "none"), 
  #                       levels=c("q1","q2","q3","q4","q5","q6","none"))
  
  q.long <- gather(q.df,group,value,
                   colnames(q.df)[4:(k+3)])
  
  
  plot <- ggplot(q.long[!is.na(q.long$site),], 
                 aes(fill=group, y=value, x=names)) + 
    geom_bar(position="stack", stat="identity", width=1) +
    facet_grid(~site, drop=TRUE, scales="free_x", shrink=TRUE, space="free",switch="x") +
    #facet_grid(~assign,drop=TRUE,space="free",scales="free_x",shrink=TRUE, switch="x") +
    #facet_grid(vars(site),drop=TRUE,space="free",scales="free_x",shrink=TRUE) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(size=6,angle=90),
          axis.ticks.x=element_blank(),
          panel.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size=18),
          legend.position="none",
          axis.text.y=element_text(size=12),
          axis.title.y=element_blank(),
          legend.text=element_text(size=16)) 
  assign(paste0("plot.k",k),plot)
  print(plot)
}

cowplot::plot_grid(plot.k2,plot.k3,
                   plot.k4,plot.k5,plot.k6,
                   plot.k7,plot.k8, ncol=1)
}
