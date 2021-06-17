#!/bin/sh

#SBATCH --account=latesgenomics
#SBATCH --time=7-00:00:00
#SBATCH --ntasks-per-node=3
#SBATCH --job-name=entropy-k%a
#SBATCH --mem=0
#SBATCH --array=1-8

module load gcc
module load gsl
module load parallel
$PATH=$PATH:/project/ltcichlidgenomics/bin/ # add path to find entropy, if necessary

if [ "$#" -ne 2]; then
    echo "two arguments passed to script; assuming these are (1) mpgl and (2) qfile"
    mpgl = $1
    qfile = $2
else if [ "$#" -ne 1 ]; then
    echo "one argument passed to script; starting entropy without a starting q"
    mpgl = $1
else
    echo "script needs either one or two arguments. exiting now."
    exit 1
fi

basename=`basename $mpgl '.entropy.mpgl'`
k=$SLURM_ARRAY_TASK_ID

export basename
export k
export mpgl
export qfile

if [ -z "$qfile" ]; then
    parallel --env basename --env k --env mpgl -delay 5 "entropy -i $mpgl -m 1 -w 0 -e 0.01 -k $k -s 10 -l 80000 -b 4000 -t 10 -o $basename.$k.{}.hdf5" ::: {1..3} # three reps, exectued in parallel
else
    parallel --env basename --env k --env qfile -delay 5 "entropy -i $mpgl -m 1 -w 0 -e 0.01 -k $k -s 10 -l 80000 -b 4000 -t 10 -q $qfile -o $basename.$k.rep{}.hdf5" ::: {1..3} # three reps, executed in parallel
fi

echo "done running all three iterations of entropy for k=$k on $basename"