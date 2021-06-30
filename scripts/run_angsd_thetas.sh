#! /bin/sh

## script for running angsd

module load gcc
module load gsl
module load htslib
module load bzip2
module load samtools

angsd_dir="/home/jrick/bin/angsd" 	# path to where angsd is installed
basedir="./" 						# directory where input files are, and where results will be written
refgenome="/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa" # reference genome used for bam alignment

bamlist=$1 #reads in bamlist name from command line

pop=`basename $bamlist .bamlist`

date
echo "beginning angsd calculation of sfs and thetas for ${pop} using reference genome ${refgenome}. all output will be written to ${basedir}angsd_output/"

## Step 1: Calculate global site allele frequency likelihoods (SAF)
# make sure to adjust parameters to those appropriate for your data
${angsd_dir}/angsd -b ${basedir}$bamlist \
	-anc $refgenome \
	-ref $refgenome \
	-doSaf 1 \
	-GL 1 \
	-doMaf 2 \
	-doMajorMinor 1 \
	-skipTriallelic 1 \
	-doGeno 9 \
	-doDepth 1 \
	-doCounts 1 \
	-doPost 1 \
	-setMinDepth 10 \
	-setMaxDepth 1000 \
	-baq 1 \
	-C 50 \
	-minMapQ 20 \
	-minQ 20 \
	-nThreads 16 \
	-out ${basedir}angsd_output/${pop}_angsd

## Step 2: Estimating max likelihood of the folded site frequency spectrum from SAF
# note: removing -nSites causes the job to be killed and maxes out memory
${angsd_dir}/misc/realSFS ${basedir}angsd_output/${pop}_angsd.saf.idx -nSites 10000000 > ${basedir}angsd_output/${pop}_angsd.sfs

## Step 2.5: Condense output to one sfs entry, using mean of windows
awk '{for(i=1; i<=NF; i++) {a[i]+=$i; if($i!="") b[i]++}}; END {for(i=1; i<=NF; i++) printf "%s%s", a[i]/b[i], (i==NF?ORS:OFS)}' ${basedir}angsd_output/${pop}_angsd.sfs > ${basedir}angsd_output/${pop}_angsd.mean.sfs

## Step 3: Estimating theta values from SFS
## if you want windowed theta values, add -win 5000 -step 1000 to do_stat call (with your desired window/step sizes)
${angsd_dir}/misc/realSFS saf2theta ${basedir}angsd_output/${pop}_angsd.saf.idx \
	-sfs ${basedir}angsd_output/${pop}_angsd.mean.sfs \
	-outname ${basedir}angsd_output/${pop}_angsd

${angsd_dir}/misc/thetaStat do_stat ${basedir}angsd_output/${pop}_angsd.thetas.idx -nChr 29 # change to your number of chromosomes

echo "finished calculating thetas for ${pop}"
date