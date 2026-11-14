#!/bin/bash


#SBATCH -J giz_param           # Job name
#SBATCH -o FDzoomb.o       # Name of stdout output file
#SBATCH -e FDzoomb.e       # Name of stderr error file
#SBATCH -p normal      # Queue (partition) name
#SBATCH -N 6               # Total # of nodes
#SBATCH -n 768              # Total # of mpi tasks
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A AST23026       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
conda deactivate

# Launch MPI code...

ibrun ./GIZMO bhrunboyuan.param
#mpirun -np 240 ./GIZMO bhrunboyuan.param        # Use ibrun instead of mpirun or mpiexec

