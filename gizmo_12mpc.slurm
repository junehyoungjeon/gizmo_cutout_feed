#!/bin/bash


#SBATCH -J giz_12edd           # Job name
#SBATCH -o FDzoom_12test.o       # Name of stdout output file
#SBATCH -e FDzoom_12test.e       # Name of stderr error file
#SBATCH -p normal      # Queue (partition) name
#SBATCH -N 5               # Total # of nodes
#SBATCH -n 280              # Total # of mpi tasks
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST22003       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
#conda deactivate

# Launch MPI code...

ibrun ./GIZMO bhrun_12mpc.param
#mpirun -np 384 ./GIZMO bhrun_12mpc.param        # Use ibrun instead of mpirun or mpiexec

