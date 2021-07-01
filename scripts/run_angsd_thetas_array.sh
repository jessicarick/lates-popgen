#! /bin/sh

## script for running angsd
## written by j. rick, june 2021
## adapted from perl script by m. haselhorst / c.a. buerkle, july 2016
## NOTE: this script has been tested for angsd version: 0.931-193-gd12895d (htslib: 1.12)
## some commands differ from previous versions of angsd! 

#SBATCH --account=latesgenomics
#SBATCH --time=3-00:00:00 		# shouldn't take more than a day, but gets slow with lots of individuals
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=124G 				# might not be necessary, depending on size of data
#SBATCH --array=0-4 			# should be one array per bamlist

module load gcc
module load gsl/2.5
module load htslib/1.12
module load bzip2/1.0.8
module load samtools/1.12

angsd_dir="/home/jrick/bin/angsd" 	# path to where angsd is installed
basedir="./" 						# directory where input files are, and where results will be written
lists=(${basedir}*.bamlist) 		# iterates through all files that end in .bamlist in the directory
refgenome="/project/latesgenomics/reference_genomes/Lates_calcarifer/lcalcarifer_genome_v3_chromosomal.fa" # reference genome used for bam alignment

bamlist=${lists[$SLURM_ARRAY_TASK_ID]}

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