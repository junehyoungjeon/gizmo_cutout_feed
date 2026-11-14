#!/bin/bash


#SBATCH -J giz_noswallow           # Job name
#SBATCH -o FDzoom.o       # Name of stdout output file
#SBATCH -e FDzoom.e       # Name of stderr error file
#SBATCH -p skx-normal      # Queue (partition) name
#SBATCH -N 5               # Total # of nodes
#SBATCH -n 240              # Total # of mpi tasks
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=junehyoungjeon@utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-astro       # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...
conda deactivate

# Launch MPI code...

mpirun -np 240 ./GIZMO bh_initial.param        # Use ibrun instead of mpirun or mpiexec

