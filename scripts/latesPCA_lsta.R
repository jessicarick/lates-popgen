###########################
# Lates stappersii only
###########################

source("packages_funcs.R")

lsta_vcfR<-read.vcfR("../data/lsta_092320_0.5_maf0.01_thin90_dp10.recode.vcf")
lsta_vcfR <- read.vcfR("../data/combined_lsta_noHets3-4_variants_0.5_maf0.01_dp5_thin90.recode.vcf")
lsta <- vcfR2genlight(lsta_vcfR)

#latesgen <- vcfR2genlight(lates_vcfR)
col.names <- unlist(strsplit(indNames(lsta),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
indNames(lsta)<-col.names.clean

#fishID <- matrix(unlist(strsplit(col.names.clean,".",fixed=TRUE)),ncol=2,byrow=T)

fishinfo <- read.csv('../../../lates_admixture/data/lates_combined_all_info_lib.csv',
                     header=TRUE,stringsAsFactors = FALSE)
colnames(fishinfo)[1] <- "New_ID"

pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                             fishinfo, by="Moran_FishID", all.x=TRUE, all.y=FALSE)

#pop(latesgen) <- pairedinfolates$entropy_dnaID

## pulling only L. stappersii
#lsta <- gl.keep.pop(latesgen_nolowcov,"Lsta")

ploidy(lsta) <- 2
lsta <- gl.compliance.check(lsta)
# lsta_filt <- gl.filter.callrate(lsta,
#                                 method="loc",
#                                 threshold=0.1,recalc=TRUE,
#                                 plot=FALSE,v=2,mono.rm=T)
lsta_filt <- lsta

lsta_info <- left_join(data.frame(Moran_FishID = indNames(lsta_filt)),fishinfo,by="Moran_FishID")
pop(lsta_filt) <- as.factor(lsta_info$Location)
lsta_alleles <- t(as.matrix(lsta_filt))
dim(lsta_alleles)
lsta_pca <- do.pca(lsta_alleles)
pcSummary_lsta <- summary(lsta_pca)
scree <- plot(lsta_pca, type="lines") 

par(mfrow=c(1,2),mar=c(5,5,1,1),oma=c(0,0,0,0))
plot(lsta_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lsta_info$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lsta_info$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lsta$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lsta$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(as.factor(lsta_info$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lsta_pca$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lsta_info$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lsta_info$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lsta$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lsta$importance[2,2]*100, 1), "%)", sep=""))
 legend("left",
        legend=levels(as.factor(lsta_info$Location)),
        col=scales::alpha(colors,0.6),
        border=NULL,pch=19,
        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)



##########################################
## removing library effects
## filtering SNPs found only in one group or the other
##########################################
pop(lsta_filt) <- as.factor(lsta_info$Library)
pop(lsta_filt)[pop(lsta_filt) == "lates02lates03"] <- "lates03"

#par(mfrow=c(1,2))

### DAPC optimization for number of axes to retain
# lsta.mat <- as.matrix(lsta_filt)
# lsta.mat.noNA <- gtools::na.replace(lsta.mat, mean, na.rm=T)
# lsta.xval <- xvalDapc(lsta.mat.noNA,
#                        grp=pop(lsta_filt),
#                        n.da=2,xval.plot=FALSE,training.set=0.9)
# lsta.xval[2:6] # number of PCs to achieve lowest MSE and highest success
# 
# dapc.lsta <- dapc(lsta_filt,
#                   n.pca = as.integer(lsta.xval[6][[1]]), 
#                   n.da = 2)
# 
# ### DAPC randomization to determine loadings threshold
# lsta.threshold <- randomize.dapc(lsta_filt,
#                                   pop=pop(lsta_filt),
#                                   npca=as.integer(lsta.xval[6][[1]]),
#                                   verbose=FALSE)
# 
# layout(matrix(c(1,2,2,2), nrow = 1, ncol = 4, byrow = TRUE))
# scatter.dapc(dapc.lsta, scree.da=FALSE, 
#              bg="white", pch=20, cell=0, 
#              cstar=0, solid=0.4, cex=2, clab=0, 
#              leg=TRUE, col=colors.vir[3:6])
# 
# loadingplot(dapc.lsta$var.contr[,1],
#             cex.lab=0.5,srt=90,byfac=FALSE,
#             xlab="SNP location", main="")
# 
# lib.loci <- locNames(lsta_filt)[dapc.lsta$var.contr > lsta.threshold]

lates01 <- gl.keep.pop(lsta_filt,
                       "lates01",mono.rm=FALSE)
lates01 <- gl.compliance.check(lates01)
# na <- glNA(lates01,alleleAsUnit = FALSE)/length(indNames(lates01))
# hist(na)
lates01.highcalls <- gl.filter.callrate(lates01,method="loc",threshold=0.5)

# lates02 <- gl.keep.pop(lsta_filt,
#                        "lates02")
# na <- glNA(lates02,alleleAsUnit = FALSE)/length(indNames(lates02))
# hist(na)
# lates02.highcalls <- lates02[,locNames(lates02)[na < 0.05]]
# 
# lates02lates03 <- gl.keep.pop(lsta_filt,
#                               "lates02lates03")
# na <- glNA(lates02lates03,alleleAsUnit = FALSE)/length(indNames(lates02lates03))
# hist(na)
# lates02lates03.highcalls <- lates02lates03[,locNames(lates02lates03)[na < 0.05]]

lates03 <- gl.keep.pop(lsta_filt,
                       "lates03", mono.rm=FALSE)
lates03@other$loc.metrics.flags <- NULL
lates03 <- gl.compliance.check(lates03)
# na <- glNA(lates03,alleleAsUnit = FALSE)/length(indNames(lates03))
# hist(na)
lates03.highcalls <- gl.filter.callrate(lates03,method="loc",threshold=0.5)

snps.all <- Reduce(intersect, list(locNames(lates01.highcalls),
                                   #locNames(lates02.highcalls),
                                   #locNames(lates02lates03.highcalls),
                                   locNames(lates03.highcalls)))

lsta.all <- gl.keep.loc(lsta_filt,snps.all)
lsta.all <- gl.compliance.check(lsta.all)

lsta.all.nolowcov <- gl.filter.callrate(lsta.all,
                                        method="ind",
                                        threshold=0.5,
                                        mono.rm=TRUE,recalc=TRUE,
                                        plot=TRUE,v=2)

# new PCA
lsta_info_nolowcov <- left_join(data.frame(Moran_FishID = indNames(lsta.all.nolowcov)),fishinfo,by="Moran_FishID")
lsta_info_nolowcov$Library[lsta_info_nolowcov$Library == "lates02lates03"] <- "lates03"

pop(lsta.all.nolowcov) <- lsta_info_nolowcov$Location
lsta_alleles_nolowcov <- t(as.matrix(lsta.all.nolowcov))
dim(lsta_alleles_nolowcov)
lsta_pca_nolowcov <- do.pca(lsta_alleles_nolowcov)
pcSummary_lsta_nolowcov <- summary(lsta_pca_nolowcov)
scree <- plot(lsta_pca_nolowcov, type="lines") 

par(mfrow=c(1,2))
plot(lsta_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[as.factor(lsta_info_nolowcov$Library)],0.5),
     bg=scales::alpha(colors[as.factor(lsta_info_nolowcov$Library)],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lsta_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lsta_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("topright",
#        legend=levels(as.factor(lsta_info_nolowcov$Library)),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(lsta_pca_nolowcov$x[,1:2],
     pch=19,
     col=scales::alpha(colors[factor(lsta_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     bg=scales::alpha(colors[factor(lsta_info_nolowcov$Location,
                                    levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),
     cex=2, xlab=paste("PC1 (", round(pcSummary_lsta_nolowcov$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_lsta_nolowcov$importance[2,2]*100, 1), "%)", sep=""))
# legend("bottom",
#        legend=levels(factor(lsta_info_nolowcov$Location,
#                             levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))),
#        col=scales::alpha(colors,0.6),
#        border=NULL,pch=19,
#        bty="n",cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

## plot distribution of missing data vs heterozygosity
lsta.missingness <- apply(t(as.matrix(lsta.all.nolowcov)),2,function(x){sum(is.na(x))/length(x)})
lsta.het <- gl.report.heterozygosity(lsta.all.nolowcov,method="ind")
par(mfrow=c(1,2))
plot(lsta.missingness,lsta.het$Ho,pch=19,
     col=scales::alpha(colors[factor(lsta_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)

plot(lsta_pca_nolowcov$x[,1],pch=19,lsta.missingness,
     ylab="Missingness",xlab="PC1 Loading",
     col=scales::alpha(colors[factor(lsta_info_nolowcov$Location,
                                     levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Cameroon","Congo","Dar"))],0.5),cex=2)
abline(v=0,lty=2,lwd=2)

##----------------------------------------------##
## IBD ANALYSES #######
#pop(latesgen_nolowcov) <- pairedinfolates$site
#lsta_gp <- genind2genpop(gl2gi(lsta))
#Dgen_lsta <- dist.genpop(lsta_gp,method=2) # Edward's distance
lsta_info_nolowcov$Location[lsta_info_nolowcov$Location == "Kipili"] <- "Kirando"
lsta_info_nolowcov$Location[lsta_info_nolowcov$Location == "Isonga"] <- "S_Mahale"

pop(lsta.all) <- lsta_info_nolowcov$Location
Dgen_lsta <- reich.fst(lsta.all,bootstrap=100,plot=TRUE,verbose=TRUE)
Dgen_lsta_long <- melt(Dgen_lsta$fsts)

locs <- read.csv("C:/Users/jrick/Dropbox/i/R_Projects/maps/tanganyika_maps/lktang_samplingLocs.csv",header=TRUE,row.names = 1)
colnames(locs) <- c("Site","Latitude","Longitude","coords")
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

lsta_locs <- locs[locs$Site %in% unique(pop(lsta.all)),]

Dgeo_lsta <- geodist::geodist(lsta_locs, paired = FALSE, sequential = FALSE, pad = FALSE,
                              measure = "geodesic")
row.names(Dgeo_lsta) <- lsta_locs$Site
colnames(Dgeo_lsta) <- lsta_locs$Site
Dgeo_lsta_long <- subset(melt(Dgeo_lsta), value!=0)

all_dist_lsta <- left_join(Dgen_lsta_long,Dgeo_lsta_long,
                           by=c("Var1","Var2"),suffix=c(".gen",".geo"))
all_dist_lsta$value.geo <- log(all_dist_lsta$value.geo / 1000)
all_dist_lsta$value.gen.std <- all_dist_lsta$value.gen/(1-all_dist_lsta$value.gen)

par(mfrow=c(1,1),xpd=FALSE)
plot(value.gen.std ~ value.geo,data=all_dist_lsta,pch=21,bg="gray",cex=2)
abline(lm(value.gen.std ~ value.geo, data=all_dist_lsta), lwd=1.5, lty=2)
summary(lm(value.gen.std ~ value.geo, data=all_dist_lsta))

## mantel test
Dgen_lsta_fst_1fst <- as.dist(pivot_wider(all_dist_lsta[,c(1,2,5)], names_from=Var2, values_from=value.gen.std)[,-c(1)], upper=TRUE, diag = TRUE)
Dgeo_lsta_km <- as.dist(pivot_wider(all_dist_lsta[,c(1,2,4)], names_from=Var2, values_from=value.geo)[,-c(1)], upper=TRUE, diag=TRUE)
ibd_lsta <- mantel.randtest(Dgen_lsta_fst_1fst,Dgeo_lsta_km)
ibd_lsta # not significant

### all plots together
colors.vir.lates <- viridis(9)
colors.rainbow.lates <- c("#195361","#33a8c7","#52e3e1","#a0e426","#fdf148","#ffab00","#f77976","#f050ae","#d883ff","#9336fd")

lsta1 <- ggplot(data=as_tibble(lsta_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
        geom_point(aes(color=as.factor(lsta_info_nolowcov$Library)),size=4, alpha=0.6) +
        xlab(paste("PC1 (", round(pcSummary_lsta_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(pcSummary_lsta_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "right",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=c("#787876","#b1b1b1"))

lsta.imp <- ggplot(data=as.tibble(t(pcSummary_lsta_nolowcov$importance)),
                   aes(y=`Proportion of Variance`,x=1:82)) +
        geom_col() + 
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank())

lsta2 <- ggplot(data=as_tibble(lsta_pca_nolowcov$x), aes(x=PC1,y=PC2)) +
                geom_point(aes(color=factor(lsta_info_nolowcov$Location,
                                            levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale","Isonga","Ikola","Mpinbwe","Kirando","Wampembe","Kasanga"))),
                           size=6, alpha=0.7) +
                xlab(paste("PC1 (", round(pcSummary_lsta_nolowcov$importance[2,1]*100, 1), "%)", sep="")) +
                ylab(paste("PC2 (", round(pcSummary_lsta_nolowcov$importance[2,2]*100, 1), "%)", sep="")) +
                theme_bw() +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title = element_text(size=18),
                      axis.text = element_text(size=14),
                      legend.position = "none",
                      legend.text = element_text(size=14),
                      legend.title = element_blank()) +
                scale_color_manual(values=colors.rainbow.lates, drop=F)  + guides(color = guide_legend(nrow = 2)) #+
        annotation_custom(ggplotGrob(lsta.imp), xmin = -0.03, xmax = 0.01, 
                          ymin = 0.05, ymax = 0.09)

lsta3 <- ggplot(data=all_dist_lsta, aes(x=value.geo, y=value.gen.std)) +
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

ggarrange(lsta2,lsta3,nrow=1)

## plot fst values as heatmap
as.data.frame(Dgen_lsta$fsts) %>% 
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
# pop(lsta.all) <- case_when(lsta_info_nolowcov$Location %in% c("Kagunga","Kigoma") ~ "North", 
#                            lsta_info_nolowcov$Location %in% c("N_Mahale","S_Mahale","Isonga") ~ "Mahale", 
#                            lsta_info_nolowcov$Location %in% c("Ikola","Mpinbwe","Kirando","Wampembe","Kasanga","Kipili") ~ "South", 
#                            TRUE ~ "Unk")
# table(pop(lsta.all))
# Dgen_lsta_basin <- reich.fst(lsta.all, bootstrap=100, plot=TRUE, verbose=TRUE)
# 
# ## plotly plot
# pcaAll.lsta <- data.frame(#sample.id = pairedinfolates$New_ID,
#         names = lsta_info_nolowcov$Moran_FishID,
#         spp = factor(lsta_info_nolowcov$entropy_dnaID),
#         fieldID = factor(lsta_info_nolowcov$Field.ID),
#         site = factor(lsta_info_nolowcov$Location,
#                       levels=c("Kagunga","Kigoma","N_Mahale","S_Mahale",
#                                "Isonga","Ikola","Mpinbwe","Kirando",
#                                "Wampembe","Kasanga","Cameroon","Congo","Dar")),
#         #country = factor(pairedinfolates$country),
#         library = factor(lsta_info_nolowcov$Library),
#         #lab = factor(pairedinfolates$lab),
#         sex = factor(lsta_info_nolowcov$Sex),
#         #missing = (pairedinfolates$Missing),
#         #heterozyg = (pairedinfolates$Heterozygosity),
#         EV1 = lsta_pca_nolowcov$x[,1],    # the first eigenvector
#         EV2 = lsta_pca_nolowcov$x[,2],    # the second eigenvector
#         EV3 = lsta_pca_nolowcov$x[,3],    # the third eigenvector
#         EV4 = lsta_pca_nolowcov$x[,4],
#         EV5 = lsta_pca_nolowcov$x[,5],
#         stringsAsFactors = FALSE)
# 
# library(plotly)
# p.lsta1 <- plot_ly(pcaAll.lsta, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir,
#                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
#         add_markers(size=I(30)) %>%
#         layout(scene = list(xaxis = list(title = 
#                                                  paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
#                                                        sep="")),
#                             yaxis = list(title = 
#                                                  paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
#                                                        sep=""))))
# p.lsta2 <- plot_ly(pcaAll.lsta, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir,
#                    text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
#         add_markers(size=I(30)) %>%
#         layout(scene = list(xaxis = list(title = 
#                                                  paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
#                                                        sep="")),
#                             yaxis = list(title = 
#                                                  paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", 
#                                                        sep=""))))
# 
# p.lsta <- subplot(p.lsta1,p.lsta2)
# p.lsta
# 
# df <- pcaAll.lsta[!is.na(pcaAll.lsta$site),]
# 
# shared_df <- SharedData$new(df)
# 
# p1 <- plot_ly(shared_df, x = ~EV1, y = ~EV2, color = ~site, colors = colors.vir, symbol = ~library,
#               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
#         add_markers(size=I(40)) %>%
#         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
#                                                        sep="")),
#                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
#                                                        sep=""))))
# p2 <- plot_ly(shared_df, x = ~EV2, y = ~EV3, color = ~site, colors = colors.vir, symbol = ~library,
#               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
#         add_markers(size=I(40)) %>%
#         layout(scene = list(xaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)",
#                                                        sep="")),
#                             yaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", 
#                                                        sep=""))))
# p3 <- plot_ly(shared_df, x = ~EV3, y = ~EV4, color = ~site, colors = colors.vir, symbol = ~library,
#               text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site), showlegend=FALSE) %>%
#         add_markers(size=I(40)) %>%
#         layout(scene = list(xaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,1]*100, 1), "%)",
#                                                        sep="")),
#                             yaxis = list(title = paste("PC4 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
#                                                        sep=""))))
# p3d <- plot_ly(shared_df, x = ~EV1, y = ~EV2, z=~EV3, color = ~site, colors = colors.lsta,  
#                text = ~paste('ID:', names, '<br>Library:', library, '<br>Site:', site)) %>%
#         add_markers() %>%
#         layout(scene = list(xaxis = list(title = paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
#                                                        sep="")),
#                             yaxis = list(title = paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
#                                                        sep="")),
#                             zaxis = list(title = paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", 
#                                                        sep=""))))
# lsta.plotly <- subplot(p1,p2,p3)
# lsta.plotly
# 
# p3d
# 
# ## beta IBD plot
# pop(lsta.all.nolowcov) <- indNames(lsta.all.nolowcov)
# Dgen_lsta_ind <- reich.fst(lsta.all.nolowcov,bootstrap=FALSE,plot=FALSE,verbose=TRUE)
# Dgen_lsta_ind_long <- melt(Dgen_lsta_ind$fsts)
# 
# locs <- read.csv("C:/Users/jrick/Dropbox/i/R_Projects/maps/tanganyika_maps/lktang_samplingLocs.csv",header=TRUE,row.names = 1)
# colnames(locs) <- c("Site","Latitude","Longitude","coords")
# locs$Site <- row.names(locs)
# locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
# locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
# locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"
# colnames(locs)[1] <- "Location"
# 
# lsta_locs_ind <- lsta_info_nolowcov %>%
#         left_join(locs, by="Location")
# 
# Dgeo_lsta_ind <- geodist::geodist(lsta_locs_ind[,c("Latitude","Longitude")], paired = FALSE, 
#                               sequential = FALSE, pad = FALSE,
#                               measure = "geodesic")
# row.names(Dgeo_lsta_ind) <- lsta_locs_ind$Moran_FishID
# colnames(Dgeo_lsta_ind) <- lsta_locs_ind$Moran_FishID
# Dgeo_lsta_ind_long <- subset(melt(Dgeo_lsta_ind), value!=0)
# 
# all_dist_lsta_ind <- left_join(Dgen_lsta_ind_long,Dgeo_lsta_ind_long,
#                            by=c("Var1","Var2"),suffix=c(".gen",".geo"))
# all_dist_lsta_ind$value.geo <- log(all_dist_lsta_ind$value.geo / 1000)
# all_dist_lsta_ind$value.gen.std <- all_dist_lsta_ind$value.gen/(1-all_dist_lsta_ind$value.gen)
# 
# plot(value.gen.std ~ value.geo, data=all_dist_lsta_ind)


