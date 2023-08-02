#PBS -l walltime=24:00:00,select=1:ncpus=1:ompthreads=1:mem=64gb
#PBS -J 1-1000

module load matlab/R2021a

cd $PBS_O_WORKDIR/../../
# matlab -batch $(basename $BASH_SOURCE .sh)
matlab -singleCompThread -batch "region_snr_primary"

exit
