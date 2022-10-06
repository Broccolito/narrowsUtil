#' @export
execute_merge_saige = function(autosomal_pattern = "FHS_EA_MRS_chrXXX.txt",
                               autosomal_only = FALSE,
                               chrx_filename = "FHS_EA_MRS_plink2_run_chrX.txt",
                               r_filename = "make_senlinplot_runner.R",
                               bash_filename = "merge_saige_runner.sh",
                               job_name = "merge_saige",
                               n_cores = 4,
                               memory_per_core = 16){

  pwd = getwd()

  while(file.exists(r_filename)){
    r_filename = gsub(pattern = ".R", replacement = "", r_filename)
    r_filename = paste0(r_filename, "_REDO.R")
  }

  while(file.exists(bash_filename)){
    bash_filename = gsub(pattern = ".sh", replacement = "", bash_filename)
    bash_filename = paste0(bash_filename, "_REDO.sh")
  }

  execution_tempalte = "module load slurm; sbatch BASH_FILE"

  bash_template = '#! /bin/bash
#SBATCH --mem-per-cpu=MEMORYG
#SBATCH --cpus-per-task=NCORES
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name JOBNAME
#SBATCH --output=JOBNAME.out
#SBATCH --partition=salem-compute

conda activate narrowsUtil
Rscript R_FILE
'

  r_template = 'library(narrowsUtil)
merge_saige_backend(autosomal_pattern = AUTOSOMAL_PATTERN,
                    autosomal_only = AUTOSOMAL_ONLY,
                    chrx_filename = CHRX_FILENAME,
                    write_file = TRUE)
'

  execution_tempalte = gsub(pattern = "BASH_FILE",
                            replacement = bash_filename,
                            execution_tempalte)

  bash_template = gsub(pattern = "R_FILE",
                       replacement = paste0(pwd, "/", r_filename),
                       bash_template) %>%
    gsub(pattern = "MEMORY", replacement = memory_per_core) %>%
    gsub(pattern = "NCORES", replacement = n_cores) %>%
    gsub(pattern = "JOBNAME", replacement = job_name)

  r_template = gsub(pattern = "AUTOSOMAL_PATTERN", replacement = autosomal_pattern, r_template) %>%
    gsub(pattern = "AUTOSOMAL_ONLY", replacement = autosomal_only) %>%
    gsub(pattern = "CHRX_FILENAME", replacement = chrx_filename)

  write(bash_template, file = bash_filename, append = FALSE)
  write(r_template, file = r_filename, append = FALSE)
  system(execution_tempalte)

}
