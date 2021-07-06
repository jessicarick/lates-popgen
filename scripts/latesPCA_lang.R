###########################
# Lates angustifrons only
# PCA and IBD analysis
#
# Script written by J. Rick, jrick@uwyo.edu
# Last updated Fall 2020
###########################

source("../packages_funcs.R")

## Import data and clean up names
lang_vcfR<-read.vcfR("../data/lang_092320_0.5_maf0.01_thin90_dp10.recode.vcf")

lang <- vcfR2genlight(lang_vcfR)
col.names <- unlist(strsplit(indNames(lang),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
indNames(lang)<-col.names.clean
ploidy(lang) <- 2
lang <- gl.compliance.check(lang)

fishinfo <- read.csv('../data/lates_combined_all_info_lib.csv',
                     header=TRUE, stringsAsFactors = FALSE)
fishinfo$Library[fishinfo$Library == "lates02lates03"] <- "lates03"

colnames(fishinfo)[1] <- "New_ID"

pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                             fishinfo, by="Moran_FishID", all.x=TRUE, all.y=FALSE)

lang_info <- left_join(data.frame(Moran_FishID = indNames(lang)),
                       fishinfo,by="Moran_FishID")

lang_info$Location[lang_info$Location == "Kipili"] <- "Kirando"
lang_info$Location[lang_info$Location == "Isonga"] <- "S_Mahale"

# Assign sampling location as population
pop(lang) <- as.factor(lang_info$Location)

#####################
# Convert to matrix and do PCA
#####################
lang_alleles <- t(as.matrix(lang))
dim(lang_alleles)
lang_pca <- do.pca(lang_alleles)
pcSummary_lang <- summary(lang_pca)
scree <- plot(lang_pca, type="lines") 

# first, plot by library
par(mfrow=c(1,2), mar=c(5,5,1,1), oma=c(0,0,0,0))
plot(lang_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lang_info$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lang_info$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lang$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lang$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lang_info$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

# now, plot by sampling location
plot(lang_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lang_info$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                      "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                      "Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lang_info$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                      "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                      "Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lang$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lang$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lang_info$Location)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)



##########################################
## removing library effects
## filtering SNPs found only in one group or the other
##########################################
pop(lang) <- as.factor(lang_info$Library)
pop(lang)[pop(lang) == "lates02lates03"] <- "lates03"
par(mfrow=c(1,2))

lates01 <- gl.keep.pop(lang,
                       "lates01", mono.rm=F)
lates01.highcalls <- gl.filter.callrate(lates01, method="loc", threshold=0.75)

lates03 <- gl.keep.pop(lang,
                       "lates03", mono.rm=F)
lates03.highcalls <- gl.filter.callrate(lates03, method="loc", threshold=0.75)

snps.all <- Reduce(intersect, list(locNames(lates01.highcalls),
                                   locNames(lates03.highcalls)))

lang.all <- gl.keep.loc(lang, snps.all)

lang.all.nolowcov <- gl.filter.callrate(lang.all,
                                        method="ind",
                                        threshold=0.5,
                                        mono.rm=TRUE, recalc=TRUE,
                                        plot=TRUE,v=2)

# new PCA
lang_info_nolowcov <- left_join(data.frame(Moran_FishID = indNames(lang.all)), 
                                          lang_info, by="Moran_FishID")
lang_info_nolowcov$Library[lang_info_nolowcov$Library == "lates02lates03"] <- "lates03"

pop(lang.all) <- lang_info_nolowcov$Location
lang_alleles_nolowcov <- t(as.matrix(lang.all))
dim(lang_alleles_nolowcov)
lang_pca_nolowcov <- do.pca(lang_alleles_nolowcov)
pcSummary_lang_nolowcov <- summary(lang_pca_nolowcov)
scree <- plot(lang_pca_nolowcov, type="lines") 

# plotted by library
par(mfrow=c(1,2), xpd=TRUE)
plot(lang_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lang_info_nolowcov$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lang_info_nolowcov$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lang_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lang_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lang_info_nolowcov$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

# plotted by sampling location
plot(lang_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lang_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                      "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                      "Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lang_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
                                      "Isonga","Ikola","Mpinbwe","Kirando","Wampembe",
                                      "Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lang_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lang_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(factor(lang_info_nolowcov$Location,
#                             levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=TRUE)

## plot distribution of missing data vs heterozygosity
lang.missingness <- apply(t(as.matrix(lang.all.nolowcov)),2,function(x){sum(is.na(x))/length(x)})
lang.het <- gl.report.heterozygosity(lang.all.nolowcov,method="ind")
par(mfrow=c(1,2))
plot(lang.missingness,lang.het$Ho,pch=19,
     col=scales::alpha(colors[factor(lang_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)

plot(lang_pca_nolowcov$x[,1],pch=19,lang.missingness,
     ylab="Missingness",xlab="PC1 Loading",
     col=scales::alpha(colors[factor(lang_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)
abline(v=0,lty=2,lwd=2)

##----------------------------------------------##
#####################
# Now, calculate genetic and geographic distances
#####################

pop(lang.all) <- lang_info_nolowcov$Location

# calculate Reich-Patterson fst between all sampling site pairs
Dgen_lang <- reich.fst(lang.all, bootstrap=100, plot=TRUE, verbose=TRUE)
Dgen_lang_long <- melt(Dgen_lang$fsts)

# import lat/long and fix site name discrepancies
locs <- read.csv("../data/lktang_samplingLocs.csv", header=TRUE, row.names = 1,
                  col.names=c("Site", "Latitude", "Longitude", "coords"))
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

# subset only those in the Lang dataset
lang_locs <- locs[locs$Site %in% unique(pop(lang.all)),]

# calcuate geographic distance between sites
Dgeo_lang <- geodist::geodist(lang_locs, paired = FALSE, sequential = FALSE, pad = FALSE,
                              measure = "geodesic")
row.names(Dgeo_lang) <- lang_locs$Site
colnames(Dgeo_lang) <- lang_locs$Site
Dgeo_lang_long <- subset(melt(Dgeo_lang), value!=0)

# Combine geographic and genetic distance, and plot
# also convert geog dist to km
# and standardize genetic dist using FST/(1-FST)
all_dist_lang <- left_join(Dgen_lang_long,Dgeo_lang_long,
                           by=c("Var1","Var2"),suffix=c(".gen",".geo"))
all_dist_lang$value.geo <- log(all_dist_lang$value.geo / 1000)
all_dist_lang$value.gen.std <- all_dist_lang$value.gen / (1-all_dist_lang$value.gen)

par(mfrow=c(1,1),xpd=FALSE)
plot(value.gen.std ~ value.geo,data=all_dist_lang,pch=21,bg="gray",cex=2)
abline(lm(value.gen.std ~ value.geo, data=all_dist_lang), lwd=1.5, lty=2)
summary(lm(value.gen.std ~ value.geo, data=all_dist_lang))

## Mantel test for IBD
Dgen_lang_fst_1fst <- as.dist(pivot_wider(all_dist_lang[,c(1,2,5)], 
                              names_from=Var2, 
                              values_from=value.gen.std)[,-c(1)], upper=TRUE)
Dgeo_lang_km <- as.dist(pivot_wider(all_dist_lang[,c(1,2,4)], 
                                    names_from=Var2, 
                                    values_from=value.geo)[,-c(1)])
ibd_lang <- mantel.randtest(Dgen_lang_fst_1fst, Dgeo_lang_km)
ibd_lang # not significant

### All plots together
colors.vir.lates <- viridis(10)
colors.rainbow.lates <- c("#195361","#33a8c7","#52e3e1",
                          "#a0e426","#fdf148","#ffab00",
                          "#f77976","#f050ae","#d883ff","#9336fd")

lang1 <- ggplot(data=as_tibble(lang_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=as.factor(lang_info_nolowcov$Library)),size=4) +
        xlab(paste("PC1 (", round(pcSummary_lang_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lang_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=c("#787876","#b1b1b1"))

lang.imp <- ggplot(data=as.tibble(t(pcSummary_lang_nolowcov$importance)),
                   aes(y=`Proportion of Variance`,x=1:ncol(pcSummary_lang_nolowcov$importance))) +
        geom_col() + 
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank())

lang2 <- ggplot(data=as_tibble(lang_pca_nolowcov$x), aes(x=V1,y=V2)) +
        geom_point(aes(color=factor(lang_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga"))),
                   size=6, alpha=0.7) +
        xlab("PC1") +
        ylab("PC2") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=colors.rainbow.lates, drop=F) #+
        annotation_custom(ggplotGrob(lang.imp), xmin = 0, xmax = 0.15, 
                          ymin = -0.17, ymax = -0.07)

lang3 <- ggplot(data=all_dist_lang, aes(x=value.geo, y=value.gen.std)) +
        geom_point(alpha=0.5, col="#787876", size=6) +
        stat_smooth(method=lm,col="#787876", fill="#787876",
                    se=TRUE,lty=2,alpha=0.2,fullrange=TRUE, size=2) +
        theme_custom() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              plot.margin = margin(30,5,10,5)) +
        ylab("Genetic Distance") +
        xlab("log Geographic Distance (km)")

ggarrange(lang2,lang3,nrow=1)

## plot fst values as heatmap
as.data.frame(Dgen_lang$fsts) %>% 
        rownames_to_column(var="pop2") %>% as_tibble() %>%
        gather(key="pop1", value="fst", -1) %>% ggplot(aes(pop1, pop2, fill= fst)) + 
        geom_tile() + 
        theme(legend.position="none") + theme_void() +
        geom_text(aes(label = round(fst, 3))) + 
        scale_fill_gradient(low="gray90",high="#900C3F") + 
        theme(axis.text = element_text(size=12), 
              axis.title = element_blank(), 
              legend.position="none") + 
        scale_x_discrete(position="top") 

## FST between basins
# pop(lang.all) <- case_when(lang_info_nolowcov$Location %in% c("Kagunga","Kigoma") ~ "North", 
#                            lang_info_nolowcov$Location %in% c("N_Mahale","S_Mahale","Isonga") ~ "Mahale", 
#                            lang_info_nolowcov$Location %in% c("Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Kipili") ~ "South", 
#                            TRUE ~ "Unk")
# table(pop(lang.all))
# Dgen_lang_basin <- reich.fst(lang.all, bootstrap=100, plot=TRUE, verbose=TRUE)