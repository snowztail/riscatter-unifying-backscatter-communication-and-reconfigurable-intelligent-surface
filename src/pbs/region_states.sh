#PBS -l walltime=72:00:00,select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -N region_states
#PBS -J 1-300

module load matlab/R2021b

cd $PBS_O_WORKDIR/..
matlab -nodesktop -nodisplay -nosplash -singleCompThread < region_states.m

exit
