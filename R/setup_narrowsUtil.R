setup_narrowsUtil = function(){

  bash_template = '#! /bin/bash
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name  narrowsUtilsetup
#SBATCH --output=narrowsUtilsetup.out
#SBATCH --partition=salem-compute

conda init bash
conda env list
conda create --name narrowsUtil --force
conda activate narrowsUtil
conda install -c conda-forge python=3.10
conda install -c conda-forge r-base
Rscript setup.R
'

  r_template = '
# install devtools if not already installed
if(!require("devtools")){
  install.packges("devtools")
  library("devtools")
}

# install narrowsUtil from GitHub
if(!require("narrowsUtil")){
  devtools::install_github("Broccolito/narrowsUtil")
  library("NarrowsUtil")
}
  '

  write(bash_template, file = "narrowsUtilsetup.sh")
  write(r_template, file = "narrowsUtilsetup.R")

  system("sbtqch narrowsUtilsetup.sh")

}
