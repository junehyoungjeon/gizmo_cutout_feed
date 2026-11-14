#!/bin/bash


#SBATCH -J gizmo           # Job name
#SBATCH -o gas.o       # Name of stdout output file
#SBATCH -e gas.e       # Name of stderr error file
#SBATCH -p skx      # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 12         # Total # of mpi tasks
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST23026       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
#conda deactivate

# Launch MPI code...
#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#ulimit -s unlimited
#ulimit -c unlimited
#export MALLOC_ARENA_MAX=1
#export I_MPI_DEBUG=5

ibrun ./GIZMO gas_sphere.param
#mpirun -np 240 ./GIZMO bhrun.param 1        # Use ibrun instead of mpirun or mpiexec

