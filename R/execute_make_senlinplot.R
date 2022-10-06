#' @export
execute_make_senlinplot = function(par_file = "parfile_path.par",
                        meta_file = "metal_output.txt",
                        snp_name = "6:160985526:G:A",
                        stats_filename = "forestplot_stats.csv",
                        senlinplot_filename = "senlinplot.png",
                        use_user_theme = FALSE,
                        user_theme = NA){

  pwd = getwd()
  r_filename = "make_senlinplot_runner.R"
  bash_filename = "make_senlinplot_runner.sh"

  if(file.exists(r_filename)){
    r_filename = gsub(pattern = ".R", replacement = "", r_filename)
    r_filename = paste0(r_filename, "_REDO.R")
  }

  if(file.exists(bash_filename)){
    bash_filename = gsub(pattern = ".sh", replacement = "", bash_filename)
    bash_filename = paste0(bash_filename, "_REDO.sh")
  }

  execution_tempalte = "module load slurm; sbatch BASH_FILE"

  bash_template = '#! /bin/bash
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -t 7-00:00
#SBATCH --job-name  make_senlin_plot
#SBATCH --output=make_senlin_plot.out
#SBATCH --partition=salem-compute

conda activate narrowsUtil
Rscript R_FILE
'

  r_template = 'library(narrowsUtil)
make_senlinplot_backend(senlinplot_stats_path = stats_filename,
                        senlinplot_filename = senlinplot_filename,
                        use_user_theme = use_user_theme,
                        user_theme = user_theme)
'

  execution_tempalte = gsub(pattern = "BASH_FILE",
                            replacement = bash_filename,
                            execution_tempalte)

  bash_template = gsub(pattern = "R_FILE",
                       replacement = paste0(pwd, "/", r_filename),
                       bash_template)

  write(bash_template, file = bash_filename, append = FALSE)
  write(r_template, file = r_filename, append = FALSE)
  system(execution_tempalte)

}
