###########################
# Lates microlepis only
###########################

source("packages_funcs.R")

lmic_vcfR<-read.vcfR("../data/lmic_092320_0.5_maf0.01_thin90_dp5.recode.vcf") 

lmic <- vcfR2genlight(lmic_vcfR)

# clean up sample names
col.names <- unlist(strsplit(indNames(lmic),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
#col.names.clean2 <- gsub(".t","",col.names.clean,fixed=TRUE)
indNames(lmic)<-col.names.clean
#fishID <- matrix(unlist(strsplit(col.names.clean,".",fixed=TRUE)),ncol=2,byrow=T)

# import sample metadata
fishinfo <- read.csv('../data/lates_combined_all_info_lib.csv',
                     header=TRUE,stringsAsFactors = FALSE)
colnames(fishinfo)[1] <- "New_ID"

pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                             fishinfo, by="Moran_FishID", all.x=TRUE, all.y=FALSE)
ploidy(lmic) <- 2
lmic <- gl.compliance.check(lmic)

lmic_info <- left_join(data.frame(Moran_FishID = indNames(lmic)),fishinfo,by="Moran_FishID")
pop(lmic) <- as.factor(lmic_info$Location)
lmic_alleles <- t(as.matrix(lmic))
dim(lmic_alleles)
lmic_pca <- do.pca(lmic_alleles)
pcSummary_lmic <- summary(lmic_pca)
scree <- plot(lmic_pca, type="lines") 

par(mfrow=c(1,4),mar=c(5,5,1,1),oma=c(0,0,0,0))
plot(lmic_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lmic_info$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lmic_info$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmic$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmic$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lmic_info$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lmic_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lmic_info$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lmic_info$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmic$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmic$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lmic_info$Location)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)


##########################################
## removing library effects
## filtering SNPs found only in one group or the other
##########################################
pop(lmic) <- as.factor(lmic_info$Library)
pop(lmic)[pop(lmic) == "lates02lates03"] <- "lates03"


lates01 <- gl.keep.pop(lmic,
                       "lates01", mono.rm=FALSE)
lates03 <- gl.keep.pop(lmic,
                       "lates03", mono.rm=FALSE)
lates01 <- gl.compliance.check(lates01)
lates03 <- gl.compliance.check(lates03)


lates01.highcalls <- gl.filter.callrate(lates01,method="loc",threshold=0.5,
                                        plot=F)
lates03.highcalls <- gl.filter.callrate(lates03,method="loc",threshold=0.5,
                                        plot=F)

snps.all <- Reduce(intersect, list(locNames(lates01.highcalls),
                                   #locNames(lates02.highcalls),
                                   #locNames(lates02lates03.highcalls),
                                   locNames(lates03.highcalls)))

lmic.all <- gl.keep.loc(lmic, snps.all)

lmic.all.nolowcov <- gl.filter.callrate(lmic.all,
                                        method="ind",
                                        threshold=0.5,
                                        mono.rm=F,recalc=FALSE,
                                        plot=TRUE,v=2)

lmic_info_nolowcov <- left_join(data.frame(Moran_FishID = indNames(lmic.all.nolowcov)),fishinfo,by="Moran_FishID")
lmic_info_nolowcov$Library[lmic_info_nolowcov$Library == "lates02lates03"] <- "lates03"


# new PCA
pop(lmic.all.nolowcov) <- lmic_info_nolowcov$Location
lmic_alleles_nolowcov <- t(as.matrix(lmic.all.nolowcov))
dim(lmic_alleles_nolowcov)
lmic_pca_nolowcov <- do.pca(lmic_alleles_nolowcov)
pcSummary_lmic_nolowcov <- summary(lmic_pca_nolowcov)
#scree <- plot(lmic_pca_nolowcov, type="lines") 

par(mfrow=c(1,2))
plot(lmic_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lmic_info_nolowcov$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lmic_info_nolowcov$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmic_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lmic_info_nolowcov$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lmic_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lmic_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lmic_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmic_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lmic_info_nolowcov$Location)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

## plot distribution of missing data vs heterozygosity
lmic.missingness <- apply(t(as.matrix(lmic.all.nolowcov)),2,function(x){sum(is.na(x))/length(x)})
lmic.het <- gl.report.heterozygosity(lmic.all.nolowcov,method="ind")
par(mfrow=c(1,2))
plot(lmic.missingness,lmic.het$Ho,pch=19,
     col=scales::alpha(colors[factor(lmic_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)

plot(lmic_pca_nolowcov$x[,1],pch=19,lmic.missingness,
     ylab="Missingness",xlab="PC1 Loading",
     col=scales::alpha(colors[factor(lmic_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)
abline(v=0,lty=2,lwd=2)

##----------------------------------------------##
#################
## IBD analyses
#################

#pop(latesgen_nolowcov) <- pairedinfolates$site
#lmic_gp <- genind2genpop(gl2gi(lmic))
#Dgen_lmic <- dist.genpop(lmic_gp,method=2) # Edward's distance

lmic_info_nolowcov$Location[lmic_info_nolowcov$Location == "Kipili"] <- "Kirando"
lmic_info_nolowcov$Location[lmic_info_nolowcov$Location == "Isonga"] <- "S_Mahale"

pop(lmic.all) <- lmic_info_nolowcov$Location
Dgen_lmic <- reich.fst(lmic.all,bootstrap=100,plot=TRUE,verbose=TRUE)
Dgen_lmic_long <- melt(Dgen_lmic$fsts)

locs <- read.csv("C:/Users/jrick/Dropbox/i/R_Projects/maps/tanganyika_maps/lktang_samplingLocs.csv",header=TRUE,row.names = 1)
colnames(locs) <- c("Site","Latitude","Longitude","coords")
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

lmic_locs <- locs[locs$Site %in% as.character(unique(pop(lmic.all))),]

Dgeo_lmic <- geodist::geodist(lmic_locs, paired = FALSE, sequential = FALSE, pad = FALSE,
                              measure = "geodesic")
row.names(Dgeo_lmic) <- lmic_locs$Site
colnames(Dgeo_lmic) <- lmic_locs$Site
Dgeo_lmic_long <- subset(melt(Dgeo_lmic), value!=0)

all_dist_lmic <- left_join(Dgen_lmic_long,Dgeo_lmic_long,
                           by=c("Var1","Var2"),suffix=c(".gen",".geo"))
all_dist_lmic$value.geo <- all_dist_lmic$value.geo / 1000
all_dist_lmic$value.gen.std <- all_dist_lmic$value.gen/(1-all_dist_lmic$value.gen)

par(mfrow=c(1,1),xpd=FALSE)
plot(value.gen.std ~ value.geo,data=all_dist_lmic,pch=21,bg="gray",cex=2)
abline(lm(value.gen.std ~ value.geo, data=all_dist_lmic), lwd=1.5, lty=2)
summary(lm(value.gen.std ~ value.geo, data=all_dist_lmic))

## mantel test
Dgen_lmic_fst_1fst <- as.dist(pivot_wider(all_dist_lmic[,c(1,2,5)], names_from=Var2, values_from=value.gen.std)[,-c(1)], upper=TRUE, diag = TRUE)
Dgeo_lmic_km <- as.dist(pivot_wider(all_dist_lmic[,c(1,2,4)], names_from=Var2, values_from=value.geo)[,-c(1)], upper=TRUE, diag=TRUE)
ibd_lmic <- mantel.randtest(Dgen_lmic_fst_1fst,Dgeo_lmic_km)
ibd_lmic # not significant

### all plots together
colors.rainbow.lates <- c("#195361","#33a8c7","#52e3e1","#a0e426","#fdf148","#ffab00","#f77976","#f050ae","#d883ff","#9336fd")

lmic1 <- ggplot(data=as_tibble(lmic_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=as.factor(lmic_info_nolowcov$Library)),size=3, alpha=0.6) +
        xlab(paste("PC1 (", round(pcSummary_lmic_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=c("#787876","#b1b1b1"))

lmic.imp <- ggplot(data=as.tibble(t(pcSummary_lmic_nolowcov$importance)),
                   aes(y=`Proportion of Variance`,x=1:33)) +
        geom_col() +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank())

lmic2 <- ggplot(data=as_tibble(lmic_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=factor(lmic_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga"))),
                   size=6, alpha=0.7) +
        xlab(paste("PC1 (", round(pcSummary_lmic_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=colors.rainbow.lates, drop=F) #+
        # annotation_custom(ggplotGrob(lmic.imp), xmin = 0.05, xmax = 0.09,
        #                   ymin = 0.035, ymax = 0.055)

lmic3 <- ggplot(data=all_dist_lmic, aes(x=value.geo, y=value.gen.std)) +
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

ggarrange(lmic2,lmic3,nrow=1)

## plot fst values as heatmap
as.data.frame(Dgen_lmic$fsts) %>% 
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

# ## FST between basins
# 
# pop(lmic.all) <- case_when(lmic_info_nolowcov$Location %in% c("Kagunga","Kigoma") ~ "North", 
#                            lmic_info_nolowcov$Location %in% c("N_Mahale","S_Mahale","Isonga") ~ "Mahale", 
#                            lmic_info_nolowcov$Location %in% c("Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Kipili") ~ "South", 
#                            TRUE ~ "Unk")
# table(pop(lmic.all))
# Dgen_lmic_basin <- reich.fst(lmic.all, bootstrap=100, plot=TRUE, verbose=TRUE)
# 
# ## plotly plot
# pcaAll.lmic <- data.frame(#sample.id = pairedinfolates$New_ID,
#         names = lmic_info_nolowcov$Moran_FishID,
#         spp = factor(lmic_info_nolowcov$entropy_dnaID),
#         fieldID = factor(lmic_info_nolowcov$Field.ID),
#         site = factor(lmic_info_nolowcov$Location,
#                       levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
#                                "Isonga","Ikola","Mpinbwe","Kirando",
#                                "Wampembe","Kasanga","Cameroon","Congo","Dar")),
#         #country = factor(pairedinfolates$country),
#         library = factor(lmic_info_nolowcov$Library),
#         #lab = factor(pairedinfolates$lab),
#         sex = factor(lmic_info_nolowcov$Sex),
#         #missing = (pairedinfolates$Missing),
#         #heterozyg = (pairedinfolates$Heterozygosity),
#         EV1 = lmic_pca_nolowcov$x[,1],    # the first eigenvector
#         EV2 = lmic_pca_nolowcov$x[,2],    # the second eigenvector
#         EV3 = lmic_pca_nolowcov$x[,3],    # the third eigenvector
#         EV4 = lmic_pca_nolowcov$x[,4],
#         EV5 = lmic_pca_nolowcov$x[,5],
#         stringsAsFactors = FALSE)
# 
# # library(plotly)
# # p.lmic1 <- plot_ly(pcaAll.lmic, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir,
# #                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title =
# #                                                  paste("PC1 (", round(pcSummary_lmic_nolowcov$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title =
# #                                                  paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p.lmic2 <- plot_ly(pcaAll.lmic, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir,
# #                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title =
# #                                                  paste("PC2 (", round(pcSummary_lmic_nolowcov$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title =
# #                                                  paste("PC3 (", round(pcSummary_lmic_nolowcov$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# #
# # p.lmic <- subplot(p.lmic1,p.lmic2)
# # p.lmic
# #
# # df <- pcaAll.lmic
# #
# # shared_df <- SharedData$new(df)
# #
# # p1 <- plot_ly(shared_df, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir, symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(50), symbol=~library) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p2 <- plot_ly(shared_df, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir, symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
# #         add_markers(size=I(50), symbol=~library) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# # p3 <- plot_ly(shared_df, x = ~EV3, y = ~EV4, color = ~site, colors = colors.vir, symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
# #         add_markers(size=I(50), symbol=~library) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC4 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p3d <- plot_ly(shared_df, x = ~EV1, y = ~EV2, z=~EV3, color = ~site, colors = colors.vir,
# #                text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers() %>%
# #         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             zaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# # lmic.plotly <- subplot(p1,p2,p3)
# # lmic.plotly
# # p3d
# #
# # ## trying some dartR functions ##
# # lmic.pcoa <- gl.pcoa(lmic.all.nolowcov)
# # gl.pcoa.plot(lmic.pcoa, lmic.all.nolowcov, labels="pop", xaxis=1, yaxis=2)
# # ggplotly()
# #
# # gl.pcoa.plot.3d(lmic.pcoa, lmic.all.nolowcov, xaxis=1, yaxis=2)
# # gl.tree.nj(lmic.all.nolowcov, type="fan")
