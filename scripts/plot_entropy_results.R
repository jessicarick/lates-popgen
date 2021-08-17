## Script for plotting entropy results for Lates popgen manuscript
## Written by Jessi Rick, last updated June 2021
## Adapted from script by Liz Mandeville

library(rhdf5)
library(abind)
library(RColorBrewer)
library(tidyverse)
library(memisc)
library(ggpubr)
source("theme_custom.R")
source("packages_funcs.R")
colors <- brewer.pal(3, "Set3")
entropy.cols <- c("#dd8d29","#e2d200","#94bf64","#46acc8","#926CA8")

## Accompanying metadata
details <- read.csv("../data/lates_all_metadata.csv",
                    header=T,na.strings = c("#N/A","n/a"))
results.dir <- "../../../results/entropy_061621/" # directory with entropy results
spp_list <- c("lang","lmar","lmic","lsta")

############################
### Part I: Plotting barplots and checking convergence
############################

for (spp in spp_list) {
  names <- read.table(paste0(results.dir,spp,".ind"), header=T,sep=",",col.names=c("ind","pop","taxon"))
  fishinfo <- left_join(names, details, by.x="ind", by.y="ind",
                        all.x=T, all.y=F)
  
  #spp <- "lsta" # used for debugging
  spp.assign <- data.frame(names=names$ind)
  
  for (i in 2:3){
    k <- i  
    
    ## Extract relevant parameters from .hdf5 format
    data1.q <- h5read(paste(results.dir,spp,"_angsd_maf0.01_dp5_miss0.5.recode.k",k,".rep1.hdf5",sep=""), "q")
    data2.q <- h5read(paste(results.dir,spp,"_angsd_maf0.01_dp5_miss0.5.recode.k",k,".rep2.hdf5",sep=""), "q")
    data3.q <- h5read(paste(results.dir,spp,"_angsd_maf0.01_dp5_miss0.5.recode.k",k,".rep3.hdf5",sep=""), "q")
    
    ## Combine chains for each parameter estimated, take mean
    allq <- abind(data1.q,data2.q,data3.q,along=1)
    q <- t(apply(allq, 2:3, mean))
    
    ## Get credible intervals and calculate mean
    q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
    q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
    print(mean(q.ci.width))
    print(median(q.ci.width))
    
    if(mean(q.ci.width) > 0.3){
      print(paste0("mean ci width is ",mean(q.ci.width)," for ",spp," at k=",k,"; only printing results for first and third rep. please check convergence."))
      allq <- abind(data1.q,data3.q,along=1)
      q <- t(apply(allq, 2:3, mean))
      
      ## Get credible intervals and calculate mean
      q.ci <- apply(allq, 2:3, quantile, probs=c(0.025,0.975))
      q.ci.width <- q.ci[2,2,]-q.ci[1,2,]
      print(mean(q.ci.width))
      print(median(q.ci.width))
      
      if(mean(q.ci.width) > 0.3) {
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
                       year=factor(fishinfo$sampling_year),
                       month=factor(fishinfo$sampling_month),
                       length=fishinfo$SL_mm,
                       juv=factor(fishinfo$juvenile),
                       #lib=factor(fishinfo$Library),
                       #plate=factor(fishinfo$Plate),
                       sex=fishinfo$pheno_sex,
                       q[,1:k]
    )
    
    q.df.max <- q.df %>%
      dplyr::select(starts_with("X")) %>%
      summarize(max = colnames(.)[max.col(.)]) %>%
      mutate(max = factor(max))
    q.df$assign <- q.df.max$max
    
    spp.assign[,k] <- q.df$assign
  
    
    # changing to long format for plotting
    q.long <- q.df %>% 
      pivot_longer(starts_with("X"),names_to="group",values_to="value")  
    
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
           # axis.text.y=element_text(size=14),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            legend.text=element_text(size=16),
            panel.border=element_blank(),
            plot.margin = unit(c(0,0.1,0,0.5), "cm")) 
    assign(paste0("plot.k",k),plot)
    print(plot)
    
    ## checking to see if assignments match with sampling year, length, juvenile, etc.
    if (k ==2) {
      q.info <- q.df %>%
        left_join(fishinfo, by=c("names" = "ind"))
      p1 <- q.info %>% 
        ggplot(aes(y=X1)) +
          geom_jitter(aes(x=X2,col=as.factor(site),shape=juvenile), size = 4, alpha=0.7, height=0.1, width=0.1) +
          theme_custom()
      p2 <- q.info %>%
        ggplot(aes(y=X1)) +
          geom_point(aes(x=SL_mm,col=as.factor(sampling_year)),size=5,alpha=0.70) +
          geom_smooth(aes(x=SL_mm),method="lm",col="black",alpha=0.2) +
          theme_custom()
      p3 <- q.info %>%
        ggplot(aes(y=X1)) +
        geom_boxplot(aes(x=as.factor(sampling_year)),alpha=0.70) +
        geom_jitter(aes(x=as.factor(sampling_year),col=assign),alpha=0.70,size=4,width=0.1,height=0.01) +
        theme_custom()
      p4 <- q.info %>%
        ggplot(aes(y=X1)) +
          geom_boxplot(aes(x=juvenile)) +
          geom_jitter(aes(x=juvenile,col=site),width=0.1,height=0.01,size=3,alpha=0.70) +
          theme_custom()
      p.all <- ggarrange(p1,p2,p3,p4,nrow=1)
      print(p.all)
    }
  }
  
  plot.all <- cowplot::plot_grid(plot.k2,plot.k3,
                     #plot.k4,plot.k5,plot.k6, 
                     ncol=1)
  plot.all2 <- ggarrange(plot.k2,plot.k3,
            #plot.k4,plot.k5,plot.k6, 
            ncol=1)
  assign(paste0("plot.all.",spp),plot.all)
  assign(paste0("plot.all2.",spp),plot.all2)
  
  assign(paste0(spp,"_assignments"),spp.assign)
  
}

ggarrange(plot.all.lsta,
          plot.all.lmic,
          plot.all.lmar,
          plot.all.lang,
          nrow=1)

ggarrange(plot.all2.lsta,
          plot.all2.lmic,
          plot.all2.lmar,
          plot.all2.lang,
          nrow=1)

entropy_assignments <- rbind(lang_assignments,
                             lmar_assignments,
                             lmic_assignments,
                             lsta_assignments)

###############################################
## Part II: results for all species combined ##
###############################################

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

## Convergence plots
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
                               levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola",
                                        "Mpinbwe","Kirando","Wampembe","Kasanga")),
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
           "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
           "#CCEBC5","#ffbc42")

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


########################################################
## Part III: Calculating intraspecific fst for groups ##
########################################################

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
#lsta_rad_gen <- vcf.rename(vcfR2genlight(lsta_vcfR))
lmic_gen <- vcf.rename(vcfR2genlight(lmic_vcfR))
lmar_gen <- vcf.rename(vcfR2genlight(lmar_vcfR))
lang_gen <- vcf.rename(vcfR2genlight(lang_vcfR))

if(sum(indNames(lang_gen) == lang_assignments$names) == length(lang_assignments$names)) {
  mean_fst_lang <- data.frame(k=seq(1:6),mean_fst=numeric(6),overlap=logical(6))
  for(k in 2:6){
    pop(lang_gen) <- lang_assignments[,k]
    fst.k <- reich.fst(lang_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
    mean_fst_lang[k,2] <- mean(fst.k$fsts,na.rm=T)
    mean_fst_lang[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
      FALSE
    } else {
      TRUE
    }
  }
} else {
  print("names do not match. try again!")
}

if(sum(indNames(lmar_gen) == lmar_assignments$names) == length(lmar_assignments$names)) {
  mean_fst_lmar <- data.frame(k=seq(1:3),mean_fst=numeric(3),overlap=logical(3))
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
  mean_fst_lmic <- data.frame(k=seq(1:3),mean_fst=numeric(3),overlap=logical(3))
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

if(sum(indNames(lsta_gen) == lsta_assignments$names) == length(lsta_assignments$names)) {
  mean_fst_lsta <- data.frame(k=seq(1:3),mean_fst=numeric(3),overlap=logical(3))
  for(k in 2:6){
    pop(lsta_gen) <- lsta_assignments[,k]
    fst.k <- reich.fst(lsta_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
    mean_fst_lsta[k,2] <- mean(fst.k$fsts,na.rm=T)
    mean_fst_lsta[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
      FALSE
    } else {
      TRUE
    }
  }
} else {
  print("names do not match. try again!")
}

if(sum(indNames(lsta_rad_gen) == lsta_rad_assignments$names) == length(lsta_rad_assignments$names)) {
  mean_fst_lsta_rad <- data.frame(k=seq(1:3),mean_fst=numeric(3),overlap=logical(3))
  for(k in 2:3){
    pop(lsta_rad_gen) <- lsta_rad_assignments[,k]
    fst.k <- reich.fst(lsta_rad_gen,bootstrap=100,plot=FALSE,verbose=TRUE)
    
    mean_fst_lsta_rad[k,2] <- mean(fst.k$fsts,na.rm=T)
    mean_fst_lsta_rad[k,3] <- if(sum(fst.k$bootstraps$min_CI < 0) < 1){
      FALSE
    } else {
      TRUE
    }
  }
} else {
  print("names do not match. try again!")
}




##############################################
### Part IV: plotting DIC for each species ###
##############################################
dic <- read_csv("../data/dic_results_072021.csv")
#colnames(dic)[1] <- "spp"

colors <- c("#2a9d8f","#e76f51","#70567d","#e9c46a","black","yellow")

min_points <- dic %>% 
  pivot_longer(col=starts_with("k"),names_to="kval",values_to="model_dic") %>%
  filter(rep == "combined_50") %>%
  group_by(spp) %>%
  slice_min(model_dic, with_ties = FALSE)


dic %>% 
  pivot_longer(col=k1:k6,names_to="kval",values_to="model_dic") %>%
  # group_by(spp,kval) %>%
  # summarize(min_dic=min(model_dic)) %>%
  # ungroup() %>%
  filter(rep %in% c("combined_50")) %>%
  #filter(spp != "Lsta-rad") %>%
  ggplot(aes(x=kval,y=model_dic/10000,col=spp,group=spp)) +
  geom_line(col="gray50",size=1.5) +
  geom_point(size=4) +
  geom_point(data=min_points,aes(x=kval,y=model_dic/10000,group=spp),size=8,col="black") +
  ylab("Model DIC (x10000)") +
  facet_wrap(~factor(spp,levels=c("Lsta","Lmic","Lmar","Lang","Lsta-rad")),scales="free",nrow=1) +
  theme_custom() +
  theme(#strip.text = element_text(size=18),
        strip.text = element_blank(),
        strip.background = element_rect(fill="white"),
        axis.text = element_text(size=18)) +
  scale_color_manual(values=colors)
