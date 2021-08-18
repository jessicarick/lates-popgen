#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=7-00:00:00
#SBATCH --ntasks-per-node=3
#SBATCH --job-name=entropy
#SBATCH --mem=0
#SBATCH --array=1-6

module load gcc
module load gsl
module load zlib
module load parallel

## script takes one or two arguments: mpgl file (required) and file with expected 
## starting values for admixture proportions (-q in entropy, optional) -- given as the basename for the files
## example usage: sbatch run_entropy_array.sh lates_all.entropy.mpgl
## example usage with qfile: sbatch run_entropy_array.sh lates_all.entropy.mpgl lates_all
## (in this example, my qfiles are lates_all_qk2.txt, lates_all_qk3.txt, etc.)

if [ "$#" -eq 2 ]; then
    echo "two arguments passed to script; assuming these are (1) mpgl and (2) qfile base"
    mpgl=$1
    qfile=$2
elif [ "$#" -eq 1 ]; then
    echo "one argument passed to script; starting entropy without a starting q"
    mpgl=$1
else
    echo "script needs either one or two arguments. exiting now."
    exit 1
fi

basename=`basename $mpgl '.mpgl'`
k=${SLURM_ARRAY_TASK_ID}

########################
## ENTROPY PARAMETERS ##
chain_length=80000 
burnin=5000 
thin=10
########################

echo "basename is $basename; k is $k"

export basename
export k
export mpgl
export qfile
export chain_length
export burnin
export thin

if [ -z "$qfile" ]; then
    parallel --line-buffer --env basename --env k --env mpgl --env chain_length --env burnin --env thin --delay 2 "echo 'starting rep {}' && /project/evolgen/bin/entropy -i $mpgl -m 1 -w 0 -e 0.01 -k $k -s 10 -l $chain_length -b $burnin -t $thin -o $basename.k$k.rep{}.hdf5" ::: {1..3} # three reps, exectued in parallel
else
    parallel --line-buffer --env basename --env k --env qfile --env chain_length --env burnin --env thin --delay 2 "echo 'starting rep {}' && /project/evolgen/bin/entropy -i $mpgl -m 1 -w 0 -e 0.01 -k $k -s 10 -l $chain_length -b $burnin -t $thin -q ${qfile}_qk${k}.txt -o $basename.k${k}.rep{}.hdf5" ::: {1..3} # three reps, executed in parallel
fi

echo "done running all three iterations of entropy for k=$k on $basename"