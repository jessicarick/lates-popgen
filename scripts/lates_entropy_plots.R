## Script for plotting entropy results for Lates popgen manuscript
## Written by Jessi Rick, last updated June 2021
## Adapted from script by Liz Mandeville

library(rhdf5)
library(abind)
library(RColorBrewer)
library(tidyverse)
library(memisc)
colors <- brewer.pal(3, "Set3")

## Accompanying metadata
details <- read.csv("../data/lates_combined_all_info_dnaID_entropy.csv",
                    header=T)
results.dir <- "../data/admixture/" # directory with entropy results

for (spp in spp_list) {
  #spp <- "lsta" # used for debugging
  for (i in 2:6){
    k <- i  
    names <- read.table(paste0(results.dir,spp,".ind"), header=F, col.names="names")
    fishinfo <- left_join(names, details, by.x="names", by.y="names",
                        all.x=T, all.y=F)
  
    ## Extract relevant parameters from .hdf5 format
    data1.q <- h5read(paste(results.dir,spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep1.hdf5",sep=""), "q")
    data2.q <- h5read(paste(results.dir,spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep2.hdf5",sep=""), "q")
    data3.q <- h5read(paste(results.dir,spp,"_092320_0.5_maf0.01_thin90.entropy.k",k,".rep3.hdf5",sep=""), "q")
  
    ## Combine chains for each parameter estimated, take mean
    allq <- abind(data1.q,data2.q,data3.q,along=1)
    q <- t(apply(allq, 2:3, mean))
  
    ## Get credible intervals and calculate mean
    q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
    q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
    print(mean(q.ci.width))
    print(median(q.ci.width))
  
    ## Plot point estimates for q
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
    q.df <- data.frame(names=fishinfo$names,
                     site=factor(as.character(fishinfo$site), 
                                 levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar")),
                     spp=factor(fishinfo$entropy_dnaID),
                     #lib=factor(fishinfo$Library),
                     #plate=factor(fishinfo$Plate),
                     #sex=fishinfo$Sex,
                     q[,1:k]
    )
  
    q.df$assign <- factor(case_when(q.df$X1 > 0.6 ~ "q1",
                           q.df$X2 > 0.6 ~ "q2",
                           #q.df$X3 > 0.6 ~ "q3",
                           #q.df$X4 > 0.6 ~ "q4",
                           # q.df$X5 > 0.6 ~ "q5",
                           TRUE ~ "none"),
                        levels=c("q1","q2",
                                 #"q3",#"q4",#"q5",#"q6",
                                 "none"))
    
    q.df <- q.df %>%
      mutate(maxq = apply(q.df,1,max)
        assign2 = colnames(q.df)[apply(q.df,1,whichmax)])

  # changing to long format for plotting
    q.long <- q.df %>% 
          pivot_longer(!c(names,site,spp),names_to="group",values_to="value")  
  
  # plotting stacked barplots
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
          strip.text.x = element_text(angle=90,size=12),
          strip.text.y = element_text(size=18),
          legend.position="none",
          axis.text.y=element_text(size=12),
          axis.title.y=element_blank(),
          legend.text=element_text(size=16)) 
    print(plot)
  }
}

##########################
## all species combined ##
##########################

k <- 4

names <- read.csv(paste0(results.dir,"lates_all.ind"), header=T)
colnames(names) <- c("names","pop","taxon")
fishinfo <- left_join(names, details, by.x="names", by.y="names",
                      all.x=T, all.y=F)

## Extract relevant parameters from .hdf5 format
data1.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90.entropy.k",k,".rep1.hdf5",sep=""), "q")
data2.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90.entropy.k",k,".rep2.hdf5",sep=""), "q")
data3.q <- h5read(paste(results.dir,"lates_all_092320_0.5_maf0.01_thin90.entropy.k",k,".rep3.hdf5",sep=""), "q")

## Combine chains for each parameter estimated, take mean
allq <- abind(data1.q,data2.q,data3.q,along=1)
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
q.df <- data.frame(names=fishinfo$names,
                   site=factor(as.character(fishinfo$site), 
                               levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar")),
                   spp=factor(fishinfo$entropy_dnaID),
                   #lib=factor(fishinfo$Library),
                   #plate=factor(fishinfo$Plate),
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

q.df <- q.df %>%
    mutate(maxq = apply(q.df,1,max)
      assign2 = colnames(q.df)[apply(q.df,1,whichmax)])

# changing to long format for plotting
q.long <- q.df %>% 
          pivot_longer(!c(names,site,spp),names_to="group",values_to="value")  

colors2<-c("#8DD3C7","#BEBADA","#FB8072","#80B1D3",
           "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#ffbc42")

plot <- ggplot(q.long[!is.na(q.long$site),], 
               aes(fill=group, y=value, x=names)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  #facet_grid(~site, drop=TRUE, scales="free_x", shrink=TRUE, space="free",switch="x") +
  facet_grid(~spp,drop=TRUE,space="free",scales="free_x",shrink=TRUE, switch="x") +
  #facet_grid(vars(site),drop=TRUE,space="free",scales="free_x",shrink=TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(angle=90,size=12),
        strip.text.y = element_text(size=18),
        legend.position="none",
        axis.text.y=element_text(size=12),
        axis.title.y=element_blank(),
        legend.text=element_text(size=16)) +
  scale_fill_manual(values=colors2[c(3,1,2,10)])
print(plot)
