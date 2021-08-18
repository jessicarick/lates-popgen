###########################
# Lates angustifrons only
# IBD analysis
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

fishinfo <- read.csv('../data/lates_all_metadata.csv',
                     header=TRUE, stringsAsFactors = FALSE)
fishinfo$Library[fishinfo$Library == "lates02lates03"] <- "lates03"

colnames(fishinfo)[1] <- "New_ID"

pairedinfolates <- left_join(data.frame(New_ID=col.names.clean),
                             fishinfo, by="New_ID", all.x=TRUE, all.y=FALSE)

lang_info <- left_join(data.frame(New_ID = indNames(lang)),
                       fishinfo,by="New_ID")

lang_info$Location[lang_info$Location == "Kipili"] <- "Kirando"
lang_info$Location[lang_info$Location == "Isonga"] <- "S_Mahale"

# Assign sampling location as population
pop(lang) <- as.factor(lang_info$Location)


#####################
# Now, calculate genetic and geographic distances
#####################

# calculate Reich-Patterson fst between all sampling site pairs
Dgen_lang <- reich.fst(lang, bootstrap=100, plot=TRUE, verbose=TRUE)
Dgen_lang_long <- melt(Dgen_lang$fsts)

# import lat/long and fix site name discrepancies
locs <- read.csv("../data/lktang_samplingLocs.csv", header=TRUE, row.names = 1,
                  col.names=c("Site", "Latitude", "Longitude", "coords"))
locs$Site <- row.names(locs)
locs$Site[locs$Site == "North Mahale"] <- "N_Mahale"
locs$Site[locs$Site == "South Mahale"] <- "S_Mahale"
locs$Site[locs$Site == "Kigoma Bay"] <- "Kigoma"

# subset only those in the Lang dataset
lang_locs <- locs[locs$Site %in% unique(pop(lang)),]

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

## import EMU PCA results
vec <- read_table2("../../../results/emu_062021/lang_092320_0.5_maf0.01_thin90_dp5.emu.eigenvecs",col_names=FALSE) # Reads in eigenvectors
val <- read_table2("../../../results/emu_062021/lang_092320_0.5_maf0.01_thin90_dp5.emu.eigenvals",col_names="eigenval") %>%
        mutate(pct = eigenval/sum(eigenval))
fam <- read_table2("../../../results/emu_062021/lang_092320_0.5_maf0.01_thin90_dp5.fam",col_names=FALSE)
colnames(fam)[1] <- "ind"

combined <- fam %>%
        dplyr::select(ind) %>%
        left_join(lang_info,by=c("ind" = "New_ID")) %>%
        bind_cols(vec) %>%
        drop_na(all_of(c("sampling_loc","final_ID")))

lang1 <- ggplot(data=combined, aes(x=X1,y=X2)) +
        geom_point(aes(color=as.factor(seq_library)),size=4) +
        xlab(paste("PC1 (", round(val$pct[1]*100, 1), "%)", sep="")) +
        ylab(paste("PC2 (", round(val$pct[2]*100, 1), "%)", sep="")) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              legend.position = "none",
              legend.text = element_text(size=14),
              legend.title = element_blank()) +
        scale_color_manual(values=c("#787876","#b1b1b1"))

lang2 <- ggplot(data=combined, aes(x=X1,y=X2)) +
        geom_point(aes(color=factor(sampling_loc,
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
              legend.title = element_blank()) 

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

## individual FSTs
pop(lang) <- lang_info$Moran_FishID
Dgen_lang_ind <- reich.fst(lang, bootstrap=FALSE, plot=FALSE, verbose=TRUE)
Dgen_lang_ind_long <- melt(Dgen_lang_ind$fsts)

