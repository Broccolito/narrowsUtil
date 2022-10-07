if(!require("dplyr")){
  install.packages("dplyr")
  library("dplyr")
}

if(!require("purrr")){
  install.packages("purrr")
  library("purrr")
}

if(!require("data.table")){
  install.packages("data.table")
  library("data.table")
}

if(!require("ggplot2")){
  install.packages("ggplot2")
  library("ggplot2")
}

if(!require("ggpubr")){
  install.packages("ggpubr")
  library("ggpubr")
}

if(!require("forestploter")){
  install.packages("forestploter")
  library("forestploter")
}

process_parfile = function(par_path = "parfile_path.par",
                           outfile_path = "metal_output.txt"){
  
  par_file = readLines(par_path)
  gwasfile = par_file[grepl(pattern = "PROCESS /", par_file)] %>%
    gsub(pattern = "PROCESS ", replacement = "")
  
  par = list()
  par[["gwasfile"]] = gwasfile
  par[["outfile"]] = outfile_path
  
  return(par)
}

get_marker_stats = function(gwas_path, snp_name = "6:160985526:G:A"){
  
  identifier_column_name = "1KG_ID"
  effectsize_column_name = "BETA"
  stderr_column_name = "SE"
  pvalue_column_name = "p.value"
  samplesize_column_name = "N"
  check_inverse_ref = TRUE
  
  cat(paste0("Read GWAS Summary Statistics from ", gwas_path, "...\n"))
  d = as.data.frame(fread(gwas_path))
  d_snp = d[d[["1KG_ID"]] == snp_name,]
  
  used_inverse = FALSE
  if(check_inverse_ref){
    if(length(unlist(strsplit(snp_name, split = ":")))!=4){
      return("SNP name not in chr:pos:ref:alt format...\n")
    }
    snp_name_inv = unlist(strsplit(snp_name, split = ":"))[c(1, 2, 4, 3)] %>%
      paste(collapse = ":")
    if(dim(d_snp)[1]==0){
      d_snp = d[d[,names(d)==identifier_column_name]==snp_name_inv,]
      used_inverse = TRUE
    }
  }
  
  d_snp = d_snp %>%
    select(all_of(identifier_column_name),
           all_of(effectsize_column_name),
           all_of(stderr_column_name),
           all_of(pvalue_column_name),
           all_of(samplesize_column_name))
  names(d_snp) = c("marker", "effect", "se", "pvalue", "n")
  
  if(used_inverse){
    cat(paste0(snp_name, " may have inversed ref and alt...\n"))
    d_snp[["effect"]] = -d_snp[["effect"]]
  }
  
  d_snp = d_snp[1,]
  d_snp = d_snp %>%
    mutate(c95_upper = effect + 1.96*se) %>%
    mutate(c95_lower = effect - 1.96*se) %>%
    mutate(pvalue = format(pvalue, scientific = TRUE, big.mark = ",", digits = 3)) %>%
    select(marker, effect, se, pvalue, c95_lower, c95_upper, n) %>%
    mutate(study = gwas_path) %>%
    select(study, everything())
  
  return(d_snp)
}

get_meta_stats = function(outfile_path, snp_name = "6:160985526:G:A"){
  
  identifier_column_name = "MarkerName"
  effectsize_column_name = "Effect"
  stderr_column_name = "StdErr"
  pvalue_column_name = "P-value"
  samplesize_column_name = "N"
  check_inverse_ref = TRUE
  
  cat(paste0("Read Meta Analysis Statistics...\n"))
  d = as.data.frame(fread(outfile_path))
  d_snp = d[d[["MarkerName"]] == snp_name,]
  
  used_inverse = FALSE
  if(check_inverse_ref){
    if(length(unlist(strsplit(snp_name, split = ":")))!=4){
      return("SNP name not in chr:pos:ref:alt format...\n")
    }
    snp_name_inv = unlist(strsplit(snp_name, split = ":"))[c(1, 2, 4, 3)] %>%
      paste(collapse = ":")
    if(dim(d_snp)[1]==0){
      d_snp = d[d[,names(d)==identifier_column_name]==snp_name_inv,]
      used_inverse = TRUE
    }
  }
  
  d_snp = d_snp %>%
    select(all_of(identifier_column_name),
           all_of(effectsize_column_name),
           all_of(stderr_column_name),
           all_of(pvalue_column_name),
           all_of(samplesize_column_name))
  names(d_snp) = c("marker", "effect", "se", "pvalue", "n")
  
  if(used_inverse){
    cat(paste0(snp_name, " may have inversed ref and alt...\n"))
    d_snp[["effect"]] = -d_snp[["effect"]]
  }
  
  d_snp = d_snp[1,]
  d_snp = d_snp %>%
    mutate(c95_upper = effect + 1.96*se) %>%
    mutate(c95_lower = effect - 1.96*se) %>%
    mutate(pvalue = format(pvalue, scientific = TRUE, big.mark = ",", digits = 3)) %>%
    select(marker, effect, se, pvalue, c95_lower, c95_upper, n) %>%
    mutate(study = "meta_analysis") %>%
    select(study, everything())
  
  return(d_snp)
}

get_senlinplot_stats = function(par_path = "parfile_path.par",
                                outfile_path = "metal_output.txt",
                                filename = "forestplot_stats.csv",
                                gwas_snp_name = "6:160985526:G:A"){
  
  par = process_parfile(par_path = par_path,
                        outfile_path = outfile_path)
  
  gwas_paths = as.list(par[["gwasfile"]])
  outfile_path = par[["outfile"]]
  gwas_stats = map(gwas_paths, get_marker_stats) %>%
    reduce(rbind.data.frame)
  meta_stats = get_meta_stats(outfile_path, snp_name = gwas_snp_name)
  
  cat(paste0("Combining GWAS Summary and Meta Analysis Statistics...\n"))
  stats = rbind.data.frame(gwas_stats, meta_stats)
  write.csv(stats, file = filename, quote = FALSE, row.names = FALSE)
  return(stats)
  
}

make_senlinplot = function(senlinplot_stats_path = "forestplot_stats.csv",
                           senlinplot_filename = "senlinplot.png"){
  
  d = read.csv(senlinplot_stats_path)
  d$study_name = strsplit(d$study, split = "[/]") %>%
    lapply(function(x){
      x[length(x)]
    }) %>%
    unlist()
  d = d %>%
    select(study_name, marker, effect, se, pvalue,
           c95_lower, c95_upper, n) %>%
    filter(!is.na(effect)) %>%
    mutate(effect = signif(effect, 3)) %>%
    mutate(se = round(se, 2)) %>%
    mutate(c95_lower = round(c95_lower, 2)) %>%
    mutate(c95_upper = round(c95_upper, 2)) %>%
    mutate(ci = paste0(effect," (", c95_lower, ", ", c95_upper, ")"))
  
  names(d) = c("Study", "marker", "Effect", "se", "P-value",
               "c95_lower","c95_upper", "N", "95% CI")
  
  d$` ` = paste(rep(" ", 20), collapse = " ")
  
  tm = forest_theme(base_size = 10,
                    refline_col = "black",
                    footnote_col = "black",
                    footnote_fontface = "bold",
                    footnote_cex = 1,
                    core = list(
                      bg_params=list(fill = rep("white",dim(d)[1])
                      )
                    )
  )
  
  plt = forest(data = select(d, Study, Effect, `P-value`, `95% CI`, N, ` `),
               est = d$Effect,
               lower = d$c95_lower,
               upper = d$c95_upper,
               ci_column = 6,
               is_summary = c(rep(FALSE, dim(d)[1]-1), TRUE),
               footnote = paste0("\n    ", d$marker[1]),
               theme = tm); plt
  
  p_wh = get_wh(plot = plt, unit = "in")
  ggsave(filename = senlinplot_filename, plot = plt,
         device = "png",
         dpi = 1200,
         width = p_wh[1], height = p_wh[2], units = "in")
  
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 4){
  par_path = as.character(args[1])
  outfile_path = as.character(args[2])
  gwas_snp_name = as.character(args[3])
  filename = as.character(args[4])
  
  stats_filename = paste0(filename, "_forestplot_stats.csv")
  plot_name = paste0(filename, "_forestplot.png")
  
  get_senlinplot_stats(par_path = par_path,
                       outfile_path = outfile_path,
                       filename = stats_filename,
                       gwas_snp_name = gwas_snp_name)
  
  make_senlinplot(senlinplot_stats_path = stats_filename,
                  senlinplot_filename = plot_name)
  
  if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
  }
  
}else{
  cat("Incorrect number of arguments supplied..\n")
}