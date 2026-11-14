#!/bin/bash


#SBATCH -J gizmo           # Job name
#SBATCH -o FDzoom_s3.o       # Name of stdout output file
#SBATCH -e FDzoom_s3.e       # Name of stderr error file
#SBATCH -p icx      # Queue (partition) name
#SBATCH -N 2               # Total # of nodes
#SBATCH -n 160              # Total # of mpi tasks
#SBATCH -t 04:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST23026       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
#conda deactivate

# Launch MPI code...
#export OMP_NUM_THREADS=1
ibrun ./GIZMO bhrunboyuan.param 1
#mpirun -np 576 ./GIZMO bhrunboyuan.param         # Use ibrun instead of mpirun or mpiexec

