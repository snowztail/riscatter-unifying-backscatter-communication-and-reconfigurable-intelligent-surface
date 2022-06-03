#PBS -l walltime=24:00:00,select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -N region_txs
#PBS -J 1-300

module load matlab/R2021a

cd $PBS_O_WORKDIR/..
matlab -nodesktop -nodisplay -nosplash -singleCompThread < region_txs.m

exit
