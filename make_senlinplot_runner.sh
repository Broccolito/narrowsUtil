#! /bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name  make_senlin_plot
#SBATCH --output=make_senlin_plot.out
#SBATCH --partition=salem-compute

conda activate narrowsUtil
Rscript /Users/wanjun/Desktop/narrowsUtil/make_senlinplot_runner.R

