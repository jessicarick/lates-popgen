###########################
# Lates mariae only
###########################

source("packages_funcs.R")

lmar_vcfR<-read.vcfR("../data/lmar_092320_0.5_maf0.01_thin90_dp10.recode.vcf") # TODO need to add to github!

lmar <- vcfR2genlight(lmar_vcfR)
col.names <- unlist(strsplit(indNames(lmar),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
#col.names.clean2 <- gsub(".t","",col.names.clean,fixed=TRUE)
indNames(lmar)<-col.names.clean
#fishID <- matrix(unlist(strsplit(col.names.clean,".",fixed=TRUE)),ncol=2,byrow=T)

fishinfo <- read.csv('../data/lates_combined_all_info_lib.csv',
                     header=TRUE,stringsAsFactors = FALSE)
colnames(fishinfo)[1] <- "New_ID"

pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                             fishinfo, by="Moran_FishID", all.x=TRUE, all.y=FALSE)
ploidy(lmar) <- 2
lmar <- gl.compliance.check(lmar)

lmar_info <- left_join(data.frame(Moran_FishID = indNames(lmar)),fishinfo,by="Moran_FishID")
pop(lmar) <- as.factor(lmar_info$Location)
lmar_alleles <- t(as.matrix(lmar))
dim(lmar_alleles)
lmar_pca <- do.pca(lmar_alleles)
pcSummary_lmar <- summary(lmar_pca)
scree <- plot(lmar_pca, type="lines") 

par(mfrow=c(1,2),mar=c(5,5,1,1),oma=c(0,0,0,0))
plot(lmar_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lmar_info$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lmar_info$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmar$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmar$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lmar_info$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lmar_pca$x[,1:2],
     pch=19, 
     col=scales::alpha(colors[factor(lmar_info$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lmar_info$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=3, xlab=paste("PC1 (", round(pcSummary_lmar$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmar$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lmar_info$Location)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

##########################################
## removing library effects
## filtering SNPs found only in one group or the other
##########################################
pop(lmar) <- as.factor(lmar_info$Library)
pop(lmar)[pop(lmar) == "lates02lates03"] <- "lates03"
#lmar@other$loc.metrics <- NULL
par(mfrow=c(1,2))

### DAPC optimization for number of axes to retain
# lmar.mat <- as.matrix(lmar)
# lmar.mat.noNA <- gtools::na.replace(lmar.mat, mean, na.rm=T)
# lmar.xval <- xvalDapc(lmar.mat.noNA,
#                        grp=pop(lmar),
#                        n.da=2,xval.plot=FALSE,training.set=0.9)
# lmar.xval[2:6] # number of PCs to achieve lowest MSE and highest success
# 
# dapc.lmar <- dapc(lmar,
#                   n.pca = as.integer(lmar.xval[6][[1]]), 
#                   n.da = 2)
# 
# ### DAPC randomization to determine loadings threshold
# lmar.threshold <- randomize.dapc(lmar,
#                                   pop=pop(lmar),
#                                   npca=as.integer(lmar.xval[6][[1]]),
#                                   verbose=FALSE)
# 
# layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
# scatter.dapc(dapc.lmar, scree.da=FALSE, 
#              bg="white", pch=20, cell=0, 
#              cstar=0, solid=0.4, cex=2, clab=0, 
#              leg=TRUE, col=colors.vir[3:6])
# 
# loadingplot(dapc.lmar$var.contr[,1],
#             cex.lab=0.5,srt=90,byfac=FALSE,
#             xlab="SNP location", main="")
# 
# lib.loci <- locNames(lmar)[dapc.lmar$var.contr > lmar.threshold]

lates01 <- gl.keep.pop(lmar,
                       "lates01", mono.rm=F)
lates01 <- gl.compliance.check(lates01)
lates03 <- gl.keep.pop(lmar,
                       "lates03", mono.rm=F)
lates03 <- gl.compliance.check(lates03)


lates01.highcalls <- gl.filter.callrate(lates01,method="loc",threshold=0.5)
lates03.highcalls <- gl.filter.callrate(lates03,method="loc",threshold=0.5)

snps.all <- Reduce(intersect, list(locNames(lates01.highcalls),
                                   #locNames(lates02.highcalls),
                                   #locNames(lates02lates03.highcalls),
                                   locNames(lates03.highcalls)))

lmar.all <- gl.keep.loc(lmar,snps.all)

lmar.all.nolowcov <- gl.filter.callrate(lmar.all,
                                        method="ind",
                                        threshold=0.5,
                                        mono.rm=TRUE,recalc=TRUE,
                                        plot=TRUE,v=2)

# new PCA
lmar_info_nolowcov <- left_join(data.frame(Moran_FishID = indNames(lmar.all.nolowcov)),fishinfo,by="Moran_FishID")
lmar_info_nolowcov$Library[lmar_info_nolowcov$Library == "lates02lates03"] <- "lates03"

pop(lmar.all.nolowcov) <- lmar_info_nolowcov$Location
lmar_alleles_nolowcov <- t(as.matrix(lmar.all.nolowcov))
dim(lmar_alleles_nolowcov)
lmar_pca_nolowcov <- do.pca(lmar_alleles_nolowcov)
pcSummary_lmar_nolowcov <- summary(lmar_pca_nolowcov)
scree <- plot(lmar_pca_nolowcov, type="lines") 

par(mfrow=c(1,2))
plot(lmar_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lmar_info_nolowcov$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lmar_info_nolowcov$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmar_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmar_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lmar_info_nolowcov$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lmar_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lmar_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lmar_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lmar_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lmar_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lmar_info_nolowcov$Location)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)


## plot distribution of missing data vs heterozygosity
lmar.missingness <- apply(t(as.matrix(lmar.all.nolowcov)),2,function(x){sum(is.na(x))/length(x)})
lmar.het <- gl.report.heterozygosity(lmar.all.nolowcov,method="ind")

par(mfrow=c(1,2))
plot(lmar.missingness,lmar.het$Ho,pch=19,
     col=scales::alpha(colors[factor(lmar_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)

plot(lmar_pca_nolowcov$x[,1],pch=19,lmar.missingness,
     ylab="Missingness",xlab="PC1 Loading",
     col=scales::alpha(colors[factor(lmar_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2)
abline(v=0,lty=2,lwd=2)

##----------------------------------------------##
################
## IBD analyses
################

#pop(latesgen_nolowcov) <- pairedinfolates$site
#lmar_gp <- genind2genpop(gl2gi(lmar))
#Dgen_lmar <- dist.genpop(lmar_gp,method=2) # Edward's distance
Dgen_lmar <- reich.fst(lmar.all.nolowcov,bootstrap=100,plot=TRUE,verbose=TRUE)
Dgen_lmar_long <- melt(Dgen_lmar$fsts)

locs <- read.csv("../data/lktang_samplingLocs.csv",header=TRUE,row.names = 1)
colnames(locs) <- c("Site","Latitude","Longitude","coords")
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

lmar_locs <- locs[locs$Site %in% unique(pop(lmar.all.nolowcov)),]

Dgeo_lmar <- geodist::geodist(lmar_locs, paired = FALSE, sequential = FALSE, pad = FALSE,
                              measure = "geodesic")
row.names(Dgeo_lmar) <- lmar_locs$Site
colnames(Dgeo_lmar) <- lmar_locs$Site
Dgeo_lmar_long <- subset(melt(Dgeo_lmar), value!=0)

all_dist_lmar <- left_join(Dgen_lmar_long,Dgeo_lmar_long,
                           by=c("Var1","Var2"),suffix=c(".gen",".geo"))
all_dist_lmar$value.geo <- all_dist_lmar$value.geo / 1000
all_dist_lmar$value.gen.std <- all_dist_lmar$value.gen/(1-all_dist_lmar$value.gen)
all_dist_lmar$value.geo.std <- log(all_dist_lmar$value.geo)

par(mfrow=c(1,1),xpd=FALSE)
plot(value.gen.std ~ value.geo.std,data=all_dist_lmar,pch=21,bg="gray",cex=2)
abline(lm(value.gen.std ~ value.geo.std, data=all_dist_lmar), lwd=1.5, lty=2)
summary(lm(value.gen.std ~ value.geo.std, data=all_dist_lmar))

## mantel test
Dgen_lmar_fst_1fst <- as.dist(pivot_wider(all_dist_lmar[,c(1,2,5)], names_from=Var2, values_from=value.gen.std)[,-c(1)], upper=TRUE, diag = TRUE)
Dgeo_lmar_km <- as.dist(pivot_wider(all_dist_lmar[,c(1,2,6)], names_from=Var2, values_from=value.geo.std)[,-c(1)], upper=TRUE, diag=TRUE)
ibd_lmar <- mantel.randtest(Dgen_lmar_fst_1fst,Dgeo_lmar_km)
ibd_lmar # not significant

### all plots together
colors.rainbow.lates <- c("#195361","#33a8c7","#52e3e1","#a0e426","#fdf148","#ffab00","#f77976","#f050ae","#d883ff","#9336fd")

lmar1 <- ggplot(data=as_tibble(lmar_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=as.factor(lmar_info_nolowcov$Library)),size=4) +
        xlab(paste("PC1 (", round(pcSummary_lmar_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lmar_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=c("#787876","#b1b1b1"))

lmar.imp <- ggplot(data=as.tibble(t(pcSummary_lmar_nolowcov$importance)),
                   aes(y=`Proportion of Variance`,x=1:ncol(pcSummary_lmar_nolowcov$importance))) +
        geom_col() +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank())

lmar2 <- ggplot(data=as_tibble(lmar_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=factor(lmar_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga"))),
                   size=6, alpha=0.7) +
        xlab(paste("PC1 (", round(pcSummary_lmar_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lmar_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=colors.rainbow.lates, drop=F) #+
        # annotation_custom(ggplotGrob(lmar.imp), xmin = -0.17, xmax = -0.07,
        #                   ymin = 0.07, ymax = 0.1)

lmar3 <- ggplot(data=all_dist_lmar, aes(x=value.geo, y=value.gen.std)) +
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

ggarrange(lmar2,lmar3,nrow=1)

## plot fst values as heatmap
as.data.frame(Dgen_lmar$fsts) %>% 
        rownames_to_column(var="pop2") %>% as_tibble() %>%
        gather(key="pop1", value="fst", -1) %>% 
        ggplot(aes(pop1, pop2, fill= fst)) + 
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
# pop(lmar.all) <- case_when(lmar_info_nolowcov$Location %in% c("Kagunga","Kigoma") ~ "North", 
#                            lmar_info_nolowcov$Location %in% c("N_Mahale","S_Mahale","Isonga") ~ "Mahale", 
#                            lmar_info_nolowcov$Location %in% c("Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Kipili") ~ "South", 
#                            TRUE ~ "Unk")
# table(pop(lmar.all))
# Dgen_lmar_basin <- reich.fst(lmar.all, bootstrap=100, plot=TRUE, verbose=TRUE)
# 
# ## plotly plot
# pcaAll.lmar <- data.frame(#sample.id = pairedinfolates$New_ID,
#         names = lmar_info_nolowcov$Moran_FishID,
#         spp = factor(lmar_info_nolowcov$entropy_dnaID),
#         fieldID = factor(lmar_info_nolowcov$Field.ID),
#         site = factor(lmar_info_nolowcov$Location,
#                       levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
#                                "Isonga","Ikola","Mpinbwe","Kirando",
#                                "Wampembe","Kasanga","Cameroon","Congo","Dar")),
#         #country = factor(pairedinfolates$country),
#         library = factor(lmar_info_nolowcov$Library),
#         #lab = factor(pairedinfolates$lab),
#         sex = factor(lmar_info_nolowcov$Sex),
#         #missing = (pairedinfolates$Missing),
#         #heterozyg = (pairedinfolates$Heterozygosity),
#         EV1 = lmar_pca_nolowcov$x[,1],    # the first eigenvector
#         EV2 = lmar_pca_nolowcov$x[,2],    # the second eigenvector
#         EV3 = lmar_pca_nolowcov$x[,3],    # the third eigenvector
#         EV4 = lmar_pca_nolowcov$x[,4],
#         EV5 = lmar_pca_nolowcov$x[,5],
#         stringsAsFactors = FALSE)
# 
# # library(plotly)
# # p.lmar1 <- plot_ly(pcaAll.lmar, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir,
# #                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title =
# #                                                  paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title =
# #                                                  paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p.lmar2 <- plot_ly(pcaAll.lmar, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir,
# #                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title =
# #                                                  paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title =
# #                                                  paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# #
# # p.lmar <- subplot(p.lmar1,p.lmar2)
# # p.lmar
# #
# # df <- pcaAll.lmar
# #
# # shared_df <- SharedData$new(df)
# #
# # p1 <- plot_ly(shared_df, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir, #symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p2 <- plot_ly(shared_df, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir, #symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# # p3 <- plot_ly(shared_df, x = ~EV3, y = ~EV4, color = ~site, colors = colors.vir, #symbol = ~library,
# #               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
# #         add_markers(size=I(30)) %>%
# #         layout(scene = list(xaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC4 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep=""))))
# # p3d <- plot_ly(shared_df, x = ~EV1, y = ~EV2, z=~EV3, color = ~site, colors = colors.lsta,
# #                text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
# #         add_markers() %>%
# #         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
# #                                                        sep="")),
# #                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
# #                                                        sep="")),
# #                             zaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)",
# #                                                        sep=""))))
# # lmar.plotly <- subplot(p1,p2,p3)
# # lmar.plotly
# # p3d
