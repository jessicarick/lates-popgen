###################
## Script for PCA and species assignment for all Lates species
##
## Analysis for Rick et al., in prep
## Written by J. Rick, jrick@uwyo.edu
## Last update: Fall 2020
###################

##################### 
## Loading packages and data
#####################
## Load required packages
source("packages_funcs.R")

## Importing the SNP data and metadata
## working with the dataset created with the following filters (from VCFtools): 
## missing data allowed = 50%, MAF = 0.01 cutoff

## Toggle for working with gbs vs rad data
type <- "gbs"

if (type == "gbs") {
  lates_vcfR <- read.vcfR("lates_all_092320_0.5_maf0.01_thin90.recode.vcf")
} else if (type == "rad") {
  lates_vcfR <- read.vcfR("combined_lsta_noHets3-4_variants_0.5_maf0.01.recode.vcf")
}

#####################
## Cleaning up the data
#####################

# convert to genlight object and define ploidy
latesgen <- vcfR2genlight(lates_vcfR)
ploidy(latesgen) <- 2
latesgen <- gl.compliance.check(latesgen)

# clean up names and import associated metadata
if (type == "gbs") {
  col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_noLates02/bamfiles/aln_"))
  col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
  col.names.clean[col.names.clean == "CEW16_135_2"] <- "CEW16_135"
  indNames(latesgen)<-col.names.clean
  
  fishinfo <- read.csv('lates_combined_all_info_lib.csv',
                     header=TRUE,stringsAsFactors = FALSE)
  
  pairedinfolates <- left_join(data.frame(Moran_FishID=col.names.clean),
                               fishinfo[,-1], by="Moran_FishID",
                               all.x=TRUE, all.y=FALSE)
  
  pairedinfolates$Library[pairedinfolates$Library == "lates02lates03"] <- "lates03"
  pop(latesgen) <- pairedinfolates$Library
  
} else if (type == "rad") {
  col.names <- unlist(strsplit(indNames(latesgen),"/project/latesgenomics/jrick/latesGBS_2018/combined_all/bwa_all_rad/aln_"))
  col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
  indNames(latesgen)<-col.names.clean
  
  fishinfo <- read.csv('lates_rad_metadata.csv',
                     header=TRUE,stringsAsFactors = FALSE)
  
  pairedinfolates <- left_join(data.frame(Moran_ID=col.names.clean),
                               fishinfo, by="Moran_ID", all.x=TRUE, all.y=FALSE)
  
  pop(latesgen) <- pairedinfolates$Library
}

# now, filter out any individuals with > 50% missing data
latesgen.nolowcov <- gl.filter.callrate(latesgen,method="ind",
                                        threshold=0.5,mono.rm=TRUE,
                                        recalc=FALSE,plot=TRUE,v=2)

#####################
## Library effects filter
## filtering SNPs found only in one group or the other
#####################

if (type == "gbs"){
  lates01 <- gl.keep.pop(latesgen.nolowcov,
                          "lates01", 
                         mono.rm = FALSE)
  lates01.highcalls <- gl.filter.callrate(lates01, 
                                          method="loc",
                                          threshold=0.75,
                                          recalc=FALSE,
                                          mono.rm = FALSE)
  lates03 <- gl.keep.pop(latesgen.nolowcov,
                          "lates03", mono.rm=FALSE)
  lates03.highcalls <- gl.filter.callrate(lates03,
                                          method="loc",
                                          threshold=0.75, 
                                          recalc=FALSE, 
                                          mono.rm=FALSE)
   
  snps.all <- Reduce(intersect, 
                     list(locNames(lates01.highcalls),
                          locNames(lates03.highcalls)))
  
} else if (type == "rad") {
  
  GQI128 <- gl.keep.pop(latesgen.nolowcov,
                          "GQI128", mono.rm = FALSE)
  GQI128.highcalls <- gl.filter.callrate(GQI128, 
                                          method="loc",
                                          threshold=0.75,
                                          recalc=FALSE,
                                          mono.rm = FALSE)

  GQI132 <- gl.keep.pop(latesgen,
                          "GQI132")
  GQI132.highcalls <- gl.filter.callrate(GQI132, 
                                          method="loc",
                                          threshold=0.75,
                                          recalc=FALSE,
                                          mono.rm = FALSE)

  GQI133 <- gl.keep.pop(latesgen.nolowcov,
                          "GQI133", mono.rm=FALSE)
  GQI133.highcalls <- gl.filter.callrate(GQI133, 
                                          method="loc",
                                          threshold=0.75,
                                          recalc=FALSE,
                                          mono.rm = FALSE)
   
  snps.all <- Reduce(intersect, list(locNames(GQI128.highcalls),
                        locNames(GQI132.highcalls),
                        locNames(GQI133.highcalls)))
}

latesgen.all <- gl.keep.loc(latesgen.nolowcov, snps.all)

# remove any indv that now have >50% missing data
latesgen_nolowcov <- gl.filter.callrate(latesgen.all,
                                          method="ind",
                                          threshold=0.5,
                                          mono.rm=T,
                                          recalc=FALSE,
                                          plot=TRUE,v=2)

## extracting genotype matrix from genlight
lates_alleles <- t(as.matrix(latesgen_nolowcov))

dim(lates_alleles)
head(colnames(lates_alleles)) ## to make sure that the names look good

#################
## Initial PCA
#################

lates_pca <- do.pca(lates_alleles)
pcSummary <- summary(lates_pca)
scree <- plot(lates_pca, type="lines") 

## Pull metadata associated with individuals in the dataset
pairedinfolates<-left_join(data.frame(Moran_FishID=indNames(latesgen_nolowcov)),
                           fishinfo,by="Moran_FishID",all.x=TRUE,all.y=FALSE) 
pairedinfolates[is.na(pairedinfolates)] <- "UNK"
head(pairedinfolates) # make sure it paired okay


## Combining fish info with PC results 

if (type == "gbs") {
  pcaAll <- data.frame(names = pairedinfolates$Moran_FishID,
                     spp = factor(pairedinfolates$entropy_dnaID),
                     fieldID = factor(pairedinfolates$Field.ID),
                     site = factor(pairedinfolates$Location),
                     library = factor(pairedinfolates$Library),
                     sex = factor(pairedinfolates$Sex),
                     EV1 = lates_pca$x[,1],    # the first eigenvector
                     EV2 = lates_pca$x[,2],    # the second eigenvector
                     EV3 = lates_pca$x[,3],    # the third eigenvector
                     EV4 = lates_pca$x[,4],
                     EV5 = lates_pca$x[,5],
                     stringsAsFactors = FALSE)
} else if (type == "rad") {
  pcaAll <- data.frame(names = pairedinfolates$Indv_ID,
                     spp = factor(pairedinfolates$Spp),
                     site = factor(pairedinfolates$Location),
                     library = factor(pairedinfolates$Library),
                     sex = factor(pairedinfolates$Sex),
                     EV1 = lates_pca$x[,1],    # the first eigenvector
                     EV2 = lates_pca$x[,2],    # the second eigenvector
                     EV3 = lates_pca$x[,3],    # the third eigenvector
                     EV4 = lates_pca$x[,4],
                     EV5 = lates_pca$x[,5],
                     stringsAsFactors = FALSE)
}

## Now, plotting the PCA results

# colored by species
par(oma=c(1,1,1,2), xpd=TRUE, mar=c(5.1, 4.1, 4.1, 1),mfrow=c(1,3))
plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))
legend("topright",legend=levels(pcaAll$spp),col=scales::alpha(colors2,0.6),border=NULL,pch=19,bty="n", cex=1, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(pcaAll$EV2, pcaAll$EV3, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))

plot(pcaAll$EV3, pcaAll$EV4, pch=21, cex=2, lwd=2, bg=scales::alpha(colors2[pcaAll$spp],0.6),col=colors2[pcaAll$spp],
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))

# colored by library (to make sure library effects are gone)
plot.new()
par(mar=c(6,6,1,4),mfrow=c(1,3),xpd=TRUE, oma=c(2,2,1,2))
plot(pcaAll$EV1, pcaAll$EV2, pch=21, cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
     col=colors[pcaAll$library],
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)
legend("bottomleft",legend=levels(pcaAll$library),col=scales::alpha(colors,0.6),border=NULL,pch=19,bty="n", cex=2, pt.cex=2, pt.lwd=2, horiz=FALSE)

plot(pcaAll$EV2, pcaAll$EV3, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5),
     col=colors[pcaAll$library],
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)

plot(pcaAll$EV3, pcaAll$EV4, pch=21,  cex=4, lwd=2, bg=scales::alpha(colors[pcaAll$library],0.5), 
     col=colors[pcaAll$library],
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
     cex.lab=3,cex.axis=2)

## 3D and interactive plot for investigating the data
library(scatterplot3d)
scatterplot3d(x=pcaAll$EV1,y=pcaAll$EV2,z=pcaAll$EV3,color=scales::alpha(c("#2a9d8f","#e76f51","#70567d","#e9c46a")[pcaAll$spp],0.4),pch=19,cex.symbols=2)

library(plotly)
p <- plot_ly(pcaAll, x = ~EV1, y = ~EV2, z = ~EV3, color = ~spp, colors = colors.vir,
             text = ~paste('ID:', names, '<br>Library:', library)) %>%
  add_markers() %>%
    
  layout(scene = list(xaxis = list(title = 
                                     paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)",
                                           sep="")),
                     yaxis = list(title = 
                                    paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", 
                                          sep="")),
                     zaxis = list(title = 
                                    paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", 
                                          sep=""))))
p

#####################
## Assign individuals to clusters
#####################

# clustering attempt using find.clusters / dapc
grp <- find.clusters(latesgen_nolowcov, max.n.clust=4)

table.value(table(pairedinfolates$entropy_dnaID, grp$grp), 
            col.lab=levels(pairedinfolates$entropy_dnaID),
            row.lab=paste("grp", 1:4))

dapc1 <- dapc(latesgen_nolowcov, grp$grp)
scatter(dapc1) # to visualize differences between groups

#####################
## Calculate Reich-Patterson FST between species
#####################

pop(latesgen_nolowcov) <- pcaAll$spp

lates_FST_spp <- reich.fst(latesgen_nolowcov,
                           bootstrap=100,
                           plot=TRUE,
                           verbose=TRUE)

# for comparison
lates_fst_spp_dartR <- gl.fst.pop(latesgen_nolowcov)
