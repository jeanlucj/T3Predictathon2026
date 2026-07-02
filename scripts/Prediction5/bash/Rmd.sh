#!/bin/bash -l
#SBATCH --time=48:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=150g
#SBATCH --tmp=150g
#SBATCH --mail-type=END
#SBATCH --mail-user=[your-email@institution.edu]
#SBATCH --output=GP.pred.RKHS.CV00
#SBATCH --job-name=GP.pred.RKHS.CV00


conda activate r_env
export RSTUDIO_PANDOC=/path/to/miniconda/envs/r_env/bin/
#export R_LIBS_USER="/path/to/R/library/"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

#R -e "rmarkdown::render('$(dirname "$0")/../analysis/combine_GRM.Rmd')"
R -e "rmarkdown::render('$(dirname "$0")/../analysis/genom_pred.Rmd')"

