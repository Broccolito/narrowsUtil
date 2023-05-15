# FILE: make_senlin_plot.R
# Created by: wagu
# Edited by m1ma
# Date: 2022/11



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

#global variable to store the effect allele from meta analysis
#meta_effectallele = "something"

#function used to process meta analysis file
get_meta_stats = function(outfile_path, snp_name = "6:160985526:G:A"){
  
  identifier_column_name = "MarkerName"
  effectsize_column_name = "Effect"
  stderr_column_name = "StdErr"
  pvalue_column_name = "P-value"
  samplesize_column_name = "N"
  check_inverse_ref = TRUE
  effectallele_column_name = "Allele1"

  cat(paste0("Read Meta Analysis Statistics...\n"))
#  if (file.exists(outfile_path)==TRUE) {
#    d = as.data.frame(fread(outfile_path))
#    d_snp = d[d[["MarkerName"]] == snp_name,]
#  }else{
#	warning(paste("File does not exist or is non-readable:", outfile_path))
#  	next
#  }
#cat("Read Meta Analysis Statistics...\n")
  if (file.exists(outfile_path)) {
  	cat("File exists:", outfile_path, "\n")
  	d <- as.data.frame(fread(outfile_path))
  	d_snp <- d[d[["MarkerName"]] == snp_name,]
  } else {
  cat(paste("File does not exist or is non-readable:", outfile_path, "\n"))
  }

  #d = as.data.frame(fread(outfile_path))
  #d_snp = d[d[["MarkerName"]] == snp_name,]
  
  #try to assigne value to effect allele
  meta_effectallele <<- d_snp[["Allele1"]]
  meta_effectallele <<- toupper(meta_effectallele)
  d_snp[["Allele1"]] = toString(meta_effectallele)
  
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
           all_of(samplesize_column_name),
	   all_of(effectallele_column_name))
  names(d_snp) = c("marker", "effect", "se", "pvalue", "n", "EA")
  
  if(used_inverse){
    cat(paste0(snp_name, " may have inversed ref and alt...\n"))
    d_snp[["effect"]] = -d_snp[["effect"]]
  }
  
  d_snp = d_snp[1,]
  d_snp = d_snp %>%
    mutate(c95_upper = effect + 1.96*se) %>%
    mutate(c95_lower = effect - 1.96*se) %>%
    mutate(pvalue = format(pvalue, scientific = TRUE, big.mark = ",", digits = 3)) %>%
    select(marker, effect, se, pvalue, c95_lower, c95_upper, n, EA) %>%
    mutate(study = "meta_analysis") %>%
    select(study, everything())

  print("meta_stats run successfully")
  return(d_snp)
}

# function used to process individual study results
get_marker_stats = function(gwas_path, snp_name = "6:160985526:G:A"){

  identifier_column_name = "1KG_ID"
  effectsize_column_name = "BETA"
  #stderr_column_name = "SE"
  #name changed bc par file has SE_gc instead of SE
  stderr_column_name = "SE_gc"
  #pvalue_column_name = "p.value"
  #name changed bc par file has P_gc insead of p.value
  pvalue_column_name = "P_gc"
  samplesize_column_name = "N"
  check_inverse_ref = TRUE

  #record effect allele
  effectallele_column_name = "EA"

  cat(paste0("Read GWAS Summary Statistics from ", gwas_path, "...\n"))
  if (file.exists(gwas_path)==TRUE) {
        cat("File exists:", gwas_path, "\n")
        d <- as.data.frame(fread(gwas_path))
        d_snp <- d[d[["1KG_ID"]] == snp_name,]


  #d = as.data.frame(fread(gwas_path))
  #d_snp = d[d[["1KG_ID"]] == snp_name,]

  #trying to find EA for individual study
  effectallele = d_snp[["EA"]]

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
           all_of(samplesize_column_name),
	   all_of(effectallele_column_name))
 
  names(d_snp) = c("marker", "effect", "se", "pvalue", "n", "EA")

  if(used_inverse){
    cat(paste0(snp_name, " may have inversed ref and alt...\n"))
    d_snp[["effect"]] = -d_snp[["effect"]]
  }

  # flip BETA value if effect allele doesn't match
  if(length(effectallele)!=0){
  	if(effectallele!=meta_effectallele){
    	d_snp[["effect"]] = d_snp[["effect"]]*(-1)
    	d_snp[["EA"]]=toString(meta_effectallele)
	}
  }


  d_snp = d_snp[1,]
  d_snp = d_snp %>%
    mutate(c95_upper = effect + 1.96*se) %>%
    mutate(c95_lower = effect - 1.96*se) %>%
    mutate(pvalue = format(pvalue, scientific = TRUE, big.mark = ",", digits = 3)) %>%
    modify_if(is.factor, as.character)%>%
    select(marker, effect, se, pvalue, c95_lower, c95_upper, n, EA) %>%
    mutate(study = gwas_path) %>%
    select(study, everything())

  rm(d)
  return(d_snp)
  } else {
  cat(paste("File does not exist or is non-readable:", gwas_path))
  }
}


get_senlinplot_stats = function(par_path = "parfile_path.par",
                                outfile_path = "metal_output.txt",
                                filename = "forestplot_stats.csv",
                                gwas_snp_name = "6:160985526:G:A"){
  par = process_parfile(par_path = par_path,
                        outfile_path = outfile_path)
  

  gwas_paths = as.list(par[["gwasfile"]])

  outfile_path = par[["outfile"]]
  meta_stats = get_meta_stats(outfile_path, snp_name = gwas_snp_name)
  gwas_stats = map(gwas_paths, get_marker_stats, snp_name = gwas_snp_name) %>%
    reduce(rbind.data.frame)
  cat(paste0("Combining GWAS Summary and Meta Analysis Statistics...\n"))
  stats = rbind.data.frame(gwas_stats, meta_stats)
  write.csv(stats, file = filename, quote = FALSE, row.names = FALSE)
  return(stats)
  
}

# function used to rename study
#rename_senlinplot_stats = function(stats_rename_file = "rename.csv",
#			original_name_file = "forestplot_stats.csv"){
#	d_rename = read.csv(stats_rename_file)
#	d_original = read.csv(original_name_file,colClasses = c("character"))
        
	#na.omit(d_original)
	#find matched study paths and new study names
	#d_rename = d_rename %>%
      		#filter(d_rename[,1] %in% d_original[,1])

	#replace the study names	
	#cat(d_original[,1])
	#cat(d_rename[,1])
#	d_original[1:(nrow(d_original)-1),1] = d_rename[,2]
#	stats = d_original %>% modify_if(is.factor, as.character)
#	write.csv(stats, file = original_name_file, quote = FALSE, row.names = FALSE)
#	return(stats)
#}

# function to rename study
rename_senlinplot_stats = function(stats_rename_file="rename.csv",
				   original_name_file= "forestplot_stats.csv"){
	d_original = read.csv(original_name_file)
	d_rename = read.csv(stats_rename_file, header=FALSE)

	# extract the first column
	original_col1 <- d_original[,1]
	rename_col1 <- d_rename[,1]

	common_rows <- intersect(original_col1,rename_col1)

	for (row in common_rows) {
		replacement <- d_rename[rename_col1 == row, 2]
		d_original[original_col1==row,1] <- replacement
	}

	write.csv(d_original,file = original_name_file,quote= FALSE, row.names=FALSE)

	cat("Study names replaced successfully. ", original_name_file, "updated.")
	#return(d_original)

}


make_senlinplot = function(senlinplot_stats_path = "forestplot_stats.csv",
                           senlinplot_filename = "senlinplot.png"){
  
  d = read.csv(senlinplot_stats_path,colClasses = c("character","character","numeric", "numeric","numeric", "numeric","numeric", "numeric","character"))
  d$study_name = strsplit(d$study, split = "[/]") %>%
    lapply(function(x){
      x[length(x)]
    }) %>%
    unlist()
  d = d %>%
    select(study_name, marker, effect, se, pvalue,
           c95_lower, c95_upper, n, EA) %>%
    filter(!is.na(effect)) %>%
    mutate(effect = signif(effect, 3)) %>%
    mutate(effect = round(effect,3)) %>%
    mutate(se = round(se, 2)) %>%
    mutate(c95_lower = round(c95_lower, 2)) %>%
    mutate(c95_upper = round(c95_upper, 2)) %>%
    mutate(ci = paste0(format(effect,drop0Trailing=F)," (",format(c95_lower,drop0Trailing=F), ", ", format(c95_upper,drop0Trailing=F), ")"))
 

  names(d) = c("Study", "marker", "Effect", "se", "P-value",
               "c95_lower","c95_upper", "N", "Effect Allele", "Effect(95% CI)")
  
  d$` ` = paste(rep(" ", 22), collapse = " ")
  
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
  
  plt = forest(data = select(d, Study, `Effect(95% CI)`, `P-value`, N, ` `),
	       est = d$Effect,
               lower = d$c95_lower,
               upper = d$c95_upper,
               ci_column = 5,
               is_summary = c(rep(FALSE, dim(d)[1]-1), TRUE),
               footnote = paste0("\n    ", d$marker[1], "; Effect Allele: ", toString(meta_effectallele)),
	       theme = tm); plt
  
  p_wh = get_wh(plot = plt, unit = "in")
  ggsave(filename = senlinplot_filename, plot = plt,
         device = "png",
         dpi = 1200,
         width = p_wh[1]+0.3, height = p_wh[2], units = "in")
  
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 5){
  par_path = as.character(args[1])
  outfile_path = as.character(args[2])
  renaming = as.character(args[3])
  gwas_snp_name = as.character(args[4])
  filename = as.character(args[5])
  
    
  stats_filename = paste0(filename, "_forestplot_stats.csv")
  plot_name = paste0(filename, "_forestplot.png")
  get_senlinplot_stats(par_path = par_path,
                       outfile_path = outfile_path,
                       filename = stats_filename,
                       gwas_snp_name = gwas_snp_name)
  
  rename_senlinplot_stats(stats_rename_file = renaming,
			  original_name_file = stats_filename)

  make_senlinplot(senlinplot_stats_path = stats_filename,
                  senlinplot_filename = plot_name)
  
  if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
  }

  rm(meta_effectallele)
  
}else{
  cat("Incorrect number of arguments supplied..\n")
  cat("Usage: Rscript make_senlin_plot_testcombine.R [Parameter file] [META ANALYSIS FILE PATH] [rename.csv] [SNPNAME] [PLOTNAME]\n")
}
