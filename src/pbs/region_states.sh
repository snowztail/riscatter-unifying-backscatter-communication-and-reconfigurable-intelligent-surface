#PBS -l walltime=72:00:00,select=1:ncpus=1:ompthreads=1:mem=4gb
#PBS -N region_states
#PBS -J 1-300
#PBS -e ${PBS_O_WORKDIR}/log/err_${PBS_JOBID}.txt
#PBS -o ${PBS_O_WORKDIR}/log/out_${PBS_JOBID}.txt

module load matlab/R2021a

cd $PBS_O_WORKDIR/..
matlab -nodesktop -nodisplay -nosplash -singleCompThread < region_states.m

exit
