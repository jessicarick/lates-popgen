## Script to extract relevant information from vcf files (PL field)
## and prepare input data for entropy in the right format 
# March 2019 -- Vivaswat Shastry
# modified June 2021 -- Jessica Rick

## USAGE: Rscript inputdataformat.R file1.vcf

library(vcfR)
library(MASS)

## flag to check if we have mixed ploidy 
mult.files<-F

args<-commandArgs(TRUE)

vcffile<-as.character(args[1])	

if(!file.exists("ploidy_inds.txt")){
    cat("Error: ploidy_inds.txt file not found! \nInclude a file titled ploidy_inds.txt with ploidy values for each ind in a new line. \nAssuming all individuals are diploid.\n")
    ploidy.file <- FALSE
} else {
    ploidy.file <- TRUE
}

if(length(args)>2){
	mult.files<-T
	ploidy1<-as.character(args[2])
	vcffile2<-as.character(args[3])
	ploidy2<-as.character(args[4])
}

## vcf to mpgl portion

print("Creating input genotype likelihood data for entropy...")

fname<-strsplit(vcffile, ".vcf")[[1]]
a.vcf<-read.vcfR(vcffile)

a.ids<-colnames(a.vcf@gt)[-1]
write.table(a.ids,file=paste0("inds_",fname,".txt"),row.names=F,col.names=F)

a.loci<-paste0(a.vcf@fix[,'CHROM'],':',a.vcf@fix[,'POS'])

a.pl<-extract.gt(a.vcf, element='PL')

## VS: maybe output some statistics from the vcf file...?
print(paste0(round(sum(is.na(a.pl))*100/(nrow(a.pl)*ncol(a.pl)),1),"% missing data in file"))
#nonNAidx<-which(apply(a.pl,1,function(x){sum(is.na(x))})==0)[1]
#
#if(is.na(nonNAidx)){
#	print("There are no loci with data for all individuals in the vcf file.")
#	print("This is OK but the ploidy file will contain a few zeros which will need to be changed.")
#	ploidy<-lapply(a.pl[10,],function(x){length(gregexpr(",",x,fixed=T)[[1]])})
#}
#else{
#	ploidy<-lapply(a.pl[nonNAidx,],function(x){length(gregexpr(",",x,fixed=T)[[1]])})
#}

if(mult.files){
	fname2<-strsplit(vcffile2, ".vcf")
	a.vcf2<-read.vcfR(vcffile2)

	a.ids2<-colnames(a.vcf2@gt)[-1]
	write.table(a.ids2,file=paste0("inds_",fname2,".txt"),row.names=F,col.names=F)

	a.pl2<-extract.gt(a.vcf2, element='PL')
	a.pl2[is.na(a.pl2)]<-paste(rep("0",ploidy2+1),collapse=",")

	a.loci2<-paste0(a.vcf2@fix[,'CHROM'],';',a.vcf2@fix[,'POS'])

	nind<-dim(a.pl)[2]+dim(a.pl2)[2]
	commonloci<-a.loci %in% a.loci2
	nloci<-dim(a.pl)[1]+dim(a.pl2)[1]-length(commonloci.idx)

	total.pl<-matrix("0", nrow=nloci, ncol=nind)
	total.pl[1:length(commonloci),1:dim(a.pl)[2]]<-a.pl[commonloci,]
	total.pl[(length(commonloci)+1):dim(a.pl)[1],1:dim(a.pl)[2]]<-a.pl[!commonloci,]

	commonloci<-a.loci2 %in% a.loci
	total.pl[1:length(commonloci),(dim(a.pl)[2]+1):nind]<-a.pl2[commonloci,]
	total.pl[(dim(a.pl)[1]+1):nloci,(dim(a.pl)[2]+1):nind]<-a.pl2[!commonloci,]

	total.pl[(length(commonloci)+1):dim(a.pl)[1],(dim(a.pl)[2]+1):nind]<-paste(rep("0",ploidy+1),collapse=",")
	total.pl[(dim(a.pl)[1]+1):nloci,1:dim(a.pl)[2]]<-paste(rep("0",ploidy2+1),collapse=",")

	total.pl<-apply(total.pl, c(1,2), function(x){as.numeric(unlist(strsplit(x,split=",")))})

	print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 

	loci.names<-a.loci[a.loci %in% a.loci2]
	loci.names<-append(a.loci[!(a.loci %in% a.loci2)],loci.names)
	loci.names<-append(a.loci2[!(a.loci2 %in% a.loci)],loci.names)

	unlink("combined.mpgl")
	cat(paste(ncol(total.pl), nrow(total.pl), "\n"), file="combined.mpgl")
	cat(paste(a.ids), file="combined.mpgl", append=T)
	for(locus in 1:nloci){
		cat(paste0("\n",loci.names[locus]," "), file="combined.mpgl", append=T)
		write.table(unlist(total.pl[locus,]), file="combined.mpgl", append=T,
					col.names=F, row.names=F, eol=" ")
	}
} else {
	a.gl<-apply(a.pl, c(1,2), function(x){as.numeric(unlist(strsplit(x,split=",")))})
	
	if (ploidy.file) {
		ploidy <- read.table("ploidy_inds.txt",header=F,stringsAsFactors=F)[,1]
	} else {
		ploidy <- 2
	}
	#ploidy<-apply(a.gl,2,function(x){length(x[[1]])-1})
    	#ploidy<-read.table("ploidy_inds.txt",header=F,stringsAsFactors=F)[,1]

	#if(length(unique(ploidy))==1){
	#	ploidy<-ploidy[1]
	#} else {
	#	print("Writing ploidy value for each individual into file...")
	#	print("[Check to make sure there are no zeros in this file]")
	#	write.table(ploidy, file="ploidy_inds.txt", row.names=F, col.names=F, sep="\n")
	#}

	nind<-dim(a.pl)[2]
	nloci<-dim(a.pl)[1]

	print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 

	for(i in 1:nrow(a.pl)){
		for(j in 1:ncol(a.pl)){
			if(sum(is.na(a.gl[,i,j]))>1){
				a.gl[,i,j]<-rep(0, ploidy[j]+1)
	        }
	    }
	}
	unlink(paste0(fname, ".mpgl"))
	cat(paste(nind, nloci, "\n"), file=paste0(fname, ".mpgl"))
	cat(paste(a.ids), file=paste0(fname, ".mpgl"), append=T)
	for(locus in 1:nloci){
		cat(paste0("\n",a.loci[locus]," "), file=paste0(fname, ".mpgl"), append=T)
		write.table(unlist(t(a.gl[,locus,])), file=paste0(fname, ".mpgl"), append=T,
					col.names=F, row.names=F, eol=" ")
	}
}


#--------------------------------------------------------------------

## mpgl to point estimates 

print("Converting genotype likelihoods into point estimates...")

if(mult.files){
	total.mean.gl<-matrix(NA, nrow=nrow(total.pl), ncol=ncol(total.pl))
	for(i in 1:nrow(total.pl)){
		for(j in 1:ncol(total.pl)){
			if(length(unique(unlist(total.pl[i,j])))!=1){
				temp<-10^(total.pl[i,j][[1]]/-10)/(sum(10^(total.pl[i,j][[1]]/-10)))
				if(length(temp)==2){
					total.mean.gl[i,j]<-sum(temp*(0:ploidy[j]))
				}
				else{
					total.mean.gl[i,j]<-sum(temp*(0:ploidy[j]))
				}
			}
		}
	}
	a.mean.gl<-total.mean.gl
}else{
	a.mean.gl<-matrix(0, nrow=nloci, ncol=nind)
	a.mean.gl<-apply(a.gl, c(2,3), function(x){temp<-sum(10^(-0.1*unlist(x)));sum((10^(-0.1*unlist(x)))*(0:max(ploidy)))/temp})
}
write.table(a.mean.gl, file="pntest_mean_gl.txt", row.names=F, col.names=F, quote=F, sep=" ")

#--------------------------------------------------------------------

# point estimates to ancestry ldak estimates

print("Obtaining admixture proportions for chain initialization...")

do.pca<-function(gmat){
	print("1. centering genotype matrix")
	gmn<-apply(gmat,1,mean, na.rm=T)
	gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat))
	gprime<-gmat-gmnmat ## remove mean

	print("2. finding pairwise covariance")
	gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
	for(i in 1:ncol(gmat)){
		for(j in i:ncol(gmat)){
			if (i==j){
				gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
			}
			else{
				gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
				gcovarmat[j,i]<-gcovarmat[i,j]
			}
		}
	}
	missing.inds<-unique(which(is.na(gcovarmat),arr.ind=T)[,2])
	print("3. final principal comp analysis")
	return(list(prcomp(x=na.omit(gcovarmat),center=TRUE,scale=FALSE), missing.inds))
}
return.val<-do.pca(a.mean.gl)
pcout<-return.val[1][[1]]
miss.inds<-return.val[2][[1]]
pcSummary<-summary(pcout)

print("4. k-means, lda and writing initial admixture proportions to file")

for(k in 2:8){
	init.admix<-matrix(0, nrow=nind, ncol=k)
	init.admix[miss.inds,]<-rep(1/k, k)

	kn<-kmeans(pcout$x[,1:10], k, iter.max=100, nstart=10, algorithm="Hartigan-Wong")
	ldakn<-lda(x=pcout$x[,1:10], grouping=kn$cluster, CV=T)
	init.admix<-ldakn$posterior
	write.table(round(init.admix, 5), quote=F, row.names=F, col.names=F, file=paste0(fname,"_qk",k,".txt"))
}

#--------------------------------------------------------------------