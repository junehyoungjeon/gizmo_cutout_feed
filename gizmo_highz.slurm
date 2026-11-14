#!/bin/bash


#SBATCH -J gizmo           # Job name
#SBATCH -o FDzoom_s3_lowrhighm.o       # Name of stdout output file
#SBATCH -e FDzoom_s3_lowrhighm.e       # Name of stderr error file
#SBATCH -p icx      # Queue (partition) name
#SBATCH -N 8               # Total # of nodes
#SBATCH -n 512              # Total # of mpi tasks
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST23026       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
#conda deactivate

# Launch MPI code...
#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#ulimit -s unlimited


ibrun ./GIZMO highz_largebox_highs8.param 
#mpirun -np 240 ./GIZMO bhrun.param 1        # Use ibrun instead of mpirun or mpiexec

