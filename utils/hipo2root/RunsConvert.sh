#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=HIPO_to_ROOT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aleksandr.bulgakov1999@gmail.com
#SBATCH --output=/farm_out/%u/%x-%A_%a-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%A_%a-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --time=24:00:00

source /u/home/manavb/myenv_clas12.sh

# Navigate to the working directory (if needed)
cd /w/hallb-scshelf2102/clas12/bulgakov/projects/momcor/proton_corr/utils/hipo2root



# Run the executable with the job array index (if needed)
srun clas12root -q -b hipo2root.c  --in=AlexSimuNewRunsEdge.dat