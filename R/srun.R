#' @export
srun = function(command, job_name = "runner",
                n_cores = 4,
                memory_per_core = 8){

  bash_template = '#! /bin/bash
#SBATCH --mem-per-cpu=MEMORYG
#SBATCH --cpus-per-task=NCORES
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name  JOBNAME
#SBATCH --output=JOBNAME.out
#SBATCH --partition=salem-compute

JOBCOMMAND
'

  bash_template = gsub(pattern = MEMORY, replacement = memory_per_core,
                       bash_template) %>%
    gsub(pattern = NCORES, replacement = n_cores) %>%
    gsub(pattern = JOBNAME, replacement = job_name) %>%
    gsub(pattern = JOBCOMMAND, replacement = command)

  write(bash_template, file = paste0(job_name, ".sh"), append = FALSE)
  system(paste0("module load slurm; sbatch ", paste0(job_name, ".sh")))

}
