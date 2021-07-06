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
  names <- read.table(paste0(results.dir,spp,".ind"), header=T,sep=",",col.names=c("ind","pop","taxon"))
  fishinfo <- left_join(names, details, by.x="ind", by.y="ind",
                        all.x=T, all.y=F)
  
  #spp <- "lsta" # used for debugging
  spp.assign <- data.frame(names=names$ind)
  
  for (i in 2:6){
    k <- i  
    
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
    
    spp.assign[,k] <- q.df$assign
    
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
            panel.border=element_blank(),
            plot.margin = unit(c(0,0.1,0,0.5), "cm")) 
    assign(paste0("plot.k",k),plot)
    print(plot)
  }
  
  plot.all <- cowplot::plot_grid(plot.k2,plot.k3,
                     plot.k4,plot.k5,plot.k6, ncol=1)
  plot.all2 <- ggarrange(plot.k2,plot.k3,
            plot.k4,plot.k5,plot.k6, ncol=1)
  assign(paste0("plot.all.",spp),plot.all)
  assign(paste0("plot.all2.",spp),plot.all2)
  
  assign(paste0(spp,"_assignments"),spp.assign)
  
}

ggarrange(plot.all.lang,
          plot.all.lmar,
          plot.all.lmic,
          #plot.all.lsta,
          nrow=1)

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

#############################################

#######
## calculating intraspecific fst for groups
########
vcf.rename <- function(x) {
  col.names <- unlist(strsplit(indNames(x),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
  col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
  col.names.clean[col.names.clean == "CEW16_135_2"] <- "CEW16_135"
  indNames(x) <- col.names.clean
  ploidy(x) <- 2
  x <- gl.compliance.check(x)
  return(x)
}

lsta_gen <- vcf.rename(vcfR2genlight(lsta_vcfR))
lmic_gen <- vcf.rename(vcfR2genlight(lmic_vcfR))
lmar_gen <- vcf.rename(vcfR2genlight(lmar_vcfR))
lang_gen <- vcf.rename(vcfR2genlight(lang_vcfR))

  if(sum(indNames(x) == lang_assignments$names) == length(lang_assignments$names)) {
    mean_fst_lang <- data.frame(k=seq(1:8),mean_fst=numeric(8),overlap=logical(8))
    for(k in 2:6){
      pop(lang_gen) <- lang_assignments[,k]
      fst.k <- reich.fst(lang_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
      mean_fst_lang[k,2] <- mean(fst.k$fsts,na.rm=T)
      mean_fst_lang[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
        FALSE
      } else {
        TRUE
      }
    } else {
      print("names do not match. try again!")
    }
  }

if(sum(indNames(lmar_gen) == lmar_assignments$names) == length(lmar_assignments$names)) {
  mean_fst_lmar <- data.frame(k=seq(1:8),mean_fst=numeric(8),overlap=logical(8))
  for(k in 2:6){
    pop(lmar_gen) <- lmar_assignments[,k]
    fst.k <- reich.fst(lmar_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
    mean_fst_lmar[k,2] <- mean(fst.k$fsts,na.rm=T)
    mean_fst_lmar[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
      FALSE
    } else {
      TRUE
    }
  }
} else {
  print("names do not match. try again!")
}

if(sum(indNames(lmic_gen) == lmic_assignments$names) == length(lmic_assignments$names)) {
  mean_fst_lmic <- data.frame(k=seq(1:8),mean_fst=numeric(8),overlap=logical(8))
  for(k in 2:6){
    pop(lmic_gen) <- lmic_assignments[,k]
    fst.k <- reich.fst(lmic_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
    mean_fst_lmic[k,2] <- mean(fst.k$fsts,na.rm=T)
    mean_fst_lmic[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
      FALSE
    } else {
      TRUE
    }
  }
} else {
  print("names do not match. try again!")
}

############################
### plotting DIC for each species
############################
dic <- read_csv("../data/dic_results_063021.csv")
colnames(dic)[1] <- "spp"

colors <- c("#2a9d8f","#70567d","#e76f51","#e9c46a")

min_points <- dic %>% 
  pivot_longer(col=starts_with("k"),names_to="kval",values_to="model_dic") %>%
  group_by(spp) %>%
  slice_min(model_dic, with_ties = FALSE)


dic %>% 
  pivot_longer(col=starts_with("k"),names_to="kval",values_to="model_dic") %>%
  ggplot(aes(x=kval,y=model_dic/10000,col=spp,group=spp)) +
  geom_line(col="gray50",size=1.5) +
  geom_point(size=4) +
  geom_point(data=min_points,aes(x=kval,y=model_dic/10000,group=spp),size=8,col="black") +
  ylab("Model DIC (x10000)") +
  facet_wrap(~factor(spp,levels=c("Lsta","Lmar","Lmic","Lang")),scales="free_y",nrow=1) +
  theme_custom() +
  theme(strip.text = element_blank()) +
  scale_color_manual(values=colors)
