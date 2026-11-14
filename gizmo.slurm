#!/bin/bash


#SBATCH -J giz_16edd4           # Job name
#SBATCH -o FDzoom.o       # Name of stdout output file
#SBATCH -e FDzoom.e       # Name of stderr error file
#SBATCH -p normal      # Queue (partition) name
#SBATCH -N 3               # Total # of nodes
#SBATCH -n 384              # Total # of mpi tasks
#SBATCH -t 04:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST23026       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
#conda deactivate

# Launch MPI code...

ibrun ./GIZMO bhrun.param 
#mpirun -np 384 ./GIZMO bhrun.param       # Use ibrun instead of mpirun or mpiexec

