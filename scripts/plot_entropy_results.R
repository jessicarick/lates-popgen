## Script for plotting entropy results for Lates popgen manuscript
## Written by Jessi Rick, last updated June 2021
## Adapted from script by Liz Mandeville

library(rhdf5)
library(abind)
library(RColorBrewer)
library(tidyverse)
library(memisc)
source("theme_custom.R")
colors <- brewer.pal(3, "Set3")

## Accompanying metadata
details <- read.csv("../data/lates_all_metadata.csv",
                    header=T)
results.dir <- "../../../lates_popgen/results/entropy_061621/" # directory with entropy results
spp_list <- c("lang","lmar","lmic")

for (spp in spp_list) {
  #spp <- "lsta" # used for debugging
  for (i in 2:6){
    k <- i  
    names <- read.table(paste0(results.dir,spp,".ind"), header=T,sep=",",col.names=c("ind","pop","taxon"))
    fishinfo <- left_join(names, details, by.x="ind", by.y="ind",
                          all.x=T, all.y=F)
    
    ## Extract relevant parameters from .hdf5 format
    data1.q <- h5read(paste(results.dir,spp,"_angsd.clean_ldak.k",k,".rep1.hdf5",sep=""), "q")[-c(1:4000),,]
    data2.q <- h5read(paste(results.dir,spp,"_angsd.clean_ldak.k",k,".rep2.hdf5",sep=""), "q")[-c(1:4000),,]
    data3.q <- h5read(paste(results.dir,spp,"_angsd.clean_ldak.k",k,".rep3.hdf5",sep=""), "q")[-c(1:4000),,]
    
    ## Combine chains for each parameter estimated, take mean
    allq <- abind(data1.q,data2.q,data3.q,along=1)
    q <- t(apply(allq, 2:3, mean))
    
    ## Get credible intervals and calculate mean
    q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
    q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
    print(mean(q.ci.width))
    print(median(q.ci.width))
    
    if(mean(q.ci.width) > 0.1){
      print(paste0("mean ci width is ",mean(q.ci.width)," for ",spp," at k=",k,"; only printing results for first and third rep. please check convergence."))
      allq <- abind(data1.q,data3.q,along=1)
      q <- t(apply(allq, 2:3, mean))
      
      ## Get credible intervals and calculate mean
      q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
      q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
      print(mean(q.ci.width))
      print(median(q.ci.width))
      
      if(mean(q.ci.width) > 0.1) {
        print(paste0("mean ci width is ",mean(q.ci.width)," for ",spp," at k=",k,"; only printing results for first rep. please check convergence."))
        allq <- abind(data1.q,along=1)
        q <- t(apply(allq, 2:3, mean))
        
        ## Get credible intervals and calculate mean
        q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
        q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
        print(mean(q.ci.width))
        print(median(q.ci.width))
      }
      #next
    }
    
    ## Plot point estimates for q -- not super helpful at k > 2
    # pdf("q_hist_lsta_k2.pdf", width=6.5, height=9)
    # par(mfrow=c(2,2))
    # for(i in 1:length(unique(fishinfo$Field.ID))){
    #   hist(q[which(fishinfo$Field.ID ==
    #                unique(fishinfo$Field.ID)[i]),1], 
    #        xlab="proportion of ancestry (q)", 
    #        main=unique(fishinfo$Field.ID)[i],
    #        col="gray90", xlim=c(0,1))
    # }
    # dev.off()
    
    ## Plot convergence plots
    pdf(paste(results.dir,spp,"_conv_plots_k",k,".pdf",sep=""))
    par(mfrow=c(4,4))
    for (i in 1:nrow(names)){
      plot(data1.q[,1,i],type="l",main=names[i,1],ylim=c(0,1))
      lines(data2.q[,1,i],col="blue")
      lines(data3.q[,1,i],col="red")
    }
    dev.off()
    
    
    # stacked barplot
    q.df <- data.frame(names=fishinfo$ind,
                       site=factor(as.character(fishinfo$sampling_loc), 
                                   levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar")),
                       spp=factor(fishinfo$final_ID),
                       #lib=factor(fishinfo$Library),
                       #plate=factor(fishinfo$Plate),
                       #sex=fishinfo$Sex,
                       q[,1:k]
    )
    
    q.df.max <- q.df %>%
      dplyr::select(starts_with("X")) %>%
      summarize(max = colnames(.)[max.col(.)]) %>%
      mutate(max = factor(max))
    q.df$assign <- q.df.max$max
    
    # q.df <- q.df %>%
    #   mutate(maxq = apply(q.df,1,max)),
    #          assign2 = colnames(q.df)[apply(q.df,1,whichmax)])
    
    # changing to long format for plotting
    q.long <- q.df %>% 
      pivot_longer(!c(names,site,spp,assign),names_to="group",values_to="value")  
    
    # plotting stacked barplots
    plot <- ggplot(q.long[!is.na(q.long$site),], 
                   aes(fill=group, y=value, x=names)) + 
      geom_bar(position="stack", stat="identity", width=1) +
      facet_grid(~site, drop=TRUE, scales="free_x", shrink=TRUE, space="free",switch="x") +
      #facet_grid(~assign,drop=TRUE,space="free",scales="free_x",shrink=TRUE, switch="x") +
      #facet_grid(vars(site),drop=TRUE,space="free",scales="free_x",shrink=TRUE) +
      theme_custom() +
      theme(axis.title.x=element_blank(),
            #axis.text.x=element_text(size=6,angle=90),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_text(size=18),
            legend.position="none",
            #axis.text.y=element_text(size=12),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            legend.text=element_text(size=16),
            panel.border=element_blank()) 
    assign(paste0("plot.k",k),plot)
    print(plot)
  }
  
  plot.all <- cowplot::plot_grid(plot.k2,plot.k3,
                     plot.k4,plot.k5,plot.k6, ncol=1)
  plot.all2 <- ggarrange(plot.k2,plot.k3,
            plot.k4,plot.k5,plot.k6, ncol=1)
  assign(paste0("plot.all.",spp),plot.all)
  assign(paste0("plot.all2.",spp),plot.all2)
  
}

##########################
## all species combined ##
##########################

k <- 4

names <- read.csv(paste0(results.dir,"lates_all.ind"), header=T)
colnames(names) <- c("ind","pop","taxon")
fishinfo <- left_join(names, details, by.x="ind", by.y="names",
                      all.x=T, all.y=F)

## Extract relevant parameters from .hdf5 format
data1.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90_dp5.k",k,".rep1.hdf5",sep=""), "q")
data2.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90_dp5.k",k,".rep2.hdf5",sep=""), "q")
data3.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90_dp5.k",k,".rep3.hdf5",sep=""), "q")

## Combine chains for each parameter estimated, take mean
allq <- abind(data1.q,data2.q,along=1)
q <- t(apply(allq, 2:3, mean))

## Get credible intervals and calculate mean
q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
print(mean(q.ci.width))
print(median(q.ci.width))

## Plot point estimates dsfor q 
# pdf("q_hist_lsta_k2.pdf", width=6.5, height=9)
# par(mfrow=c(2,2))
# for(i in 1:length(unique(fishinfo$Field.ID))){
#   hist(q[which(fishinfo$Field.ID ==
#                unique(fishinfo$Field.ID)[i]),1], 
#        xlab="proportion of ancestry (q)", 
#        main=unique(fishinfo$Field.ID)[i],
#        col="gray90", xlim=c(0,1))
# }
# dev.off()

pdf(paste(results.dir,"lates_all_conv_plots_k",k,".pdf",sep=""))
par(mfrow=c(4,4))
for (i in 1:nrow(names)){
  plot(data1.q[,1,i],type="l",main=names[i,1],ylim=c(0,1))
  lines(data2.q[,1,i],col="blue")
  lines(data3.q[,1,i],col="red")
}
dev.off()


# stacked barplot
q.df <- data.frame(names=fishinfo$ind,
                   site=factor(as.character(fishinfo$sampling_loc), 
                               levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar")),
                   spp=factor(fishinfo$final_ID),
                   #lib=factor(fishinfo$Library),
                   #plate=factor(fishinfo$Plate),
                   #sex=fishinfo$Sex,
                   q[,1:k]
)

q.df.max <- q.df %>%
  dplyr::select(starts_with("X")) %>%
  summarize(max = colnames(.)[max.col(.)]) %>%
  mutate(max = factor(max))
q.df$assign <- q.df.max$max

# changing to long format for plotting
q.long <- q.df %>% 
  pivot_longer(!c(names,site,spp),names_to="group",values_to="value")  

colors2<-c("#8DD3C7","#BEBADA","#FB8072","#80B1D3",
           "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#ffbc42")

anno_lines <- q.long %>%
  group_by(spp) %>%
  drop_na() %>%
  tally() %>%
  mutate(xmin=1,
            xmax=n,
            ymin=1.03,
            ymax=1.03)

plot <- ggplot(q.long[!is.na(q.long$site),], 
               aes(fill=group, y=value, x=names)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  #facet_grid(~site, drop=TRUE, scales="free_x", shrink=TRUE, space="free",switch="x") +
  facet_grid(~spp,drop=TRUE,space="free",scales="free_x",shrink=TRUE, switch="x") +
  #facet_grid(vars(site),drop=TRUE,space="free",scales="free_x",shrink=TRUE) +
  theme_custom() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        #strip.text.x = element_blank(),
        strip.text.y = element_text(size=18),
        legend.position="none",
        axis.text.y=element_text(size=12),
        axis.title.y=element_blank(),
        legend.text=element_text(size=16),
        panel.border = element_blank()) +
  scale_fill_manual(values=colors2[c(10,1,3,2)]) +
  geom_segment(data=anno_lines,mapping=aes(x=xmin-0.5,xend=xmax/4+0.5,y=ymin,yend=ymax),lwd=3,inherit.aes=FALSE)
print(plot)
