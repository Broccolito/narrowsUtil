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

if(!require("knitr")){
  install.packages("knitr")
  library("knitr")
}

merge_saige = function(autosomal_pattern = "FHS_EA_MRS_chrXXX.txt",
                       autosomal_only = FALSE,
                       chrx_filename = "FHS_EA_MRS_plink2_run_chrX.txt",
                       write_file = TRUE){

  files = as.list(
    sapply(1:22, function(x){
      gsub("XXX", x, autosomal_pattern)
    })
  )

  if(autosomal_only){
    all_files = c(unlist(files))
  }else{
    all_files = c(unlist(files), chrx_filename)
  }

  sanity_check_table = tibble(
    File = all_files,
    Found = all_files %in% list.files()
  )
  print(kable(sanity_check_table))
  if(any(!sanity_check_table$Found)){
    return("Some of the GWAS results files specified do not exist...\n")
  }else{
    cat("\nAll GWAS results files are found...\n\n")
  }

  filter_autosomal_tables = function(f){
    d = as.data.frame(fread(f))
    d_filtered = d %>%
      filter(imputationInfo > 0.3) %>%
      filter(AC_Allele2 > 10) %>%
      filter(SE > 0) %>%
      filter(p.value >= 0, p.value <= 1) %>%
      filter(BETA != 0)
    cat(paste0(f, " Merged...\n"))
    return(d_filtered)
  }

  filter_chrx_table = function(chrx_filename, sample_size){
    d = as.data.frame(fread(chrx_filename))
    chrx_table = tibble(CHR = 23, POS = d$POS, SNPID = d$ID,
                        Allele1 = d$REF, Allele2 = d$ALT,
                        AF_Allele2 = d$A1_FREQ, N = sample_size,
                        BETA = d$BETA, SE = d$SE, Tstat = d$T_STAT, p.value = d$P) %>%
      mutate(AC_Allele2 = sample_size*AF_Allele2) %>%
      filter(AC_Allele2 > 10) %>%
      filter(SE > 0) %>%
      filter(p.value >= 0, p.value <= 1) %>%
      filter(BETA != 0)
    cat(paste0(chrx_filename, " Merged...\n"))
    return(chrx_table)
  }

  # Filter autosomal table
  autosomal_table = map(files, filter_autosomal_tables)
  autosomal_table = autosomal_table %>%
    reduce(rbind.data.frame)

  if(!autosomal_only){
    # Filter sex chromosome table
    sample_size = max(autosomal_table$N)
    chrx_table = filter_chrx_table(chrx_filename, sample_size)
    chrx_table = full_join(head(autosomal_table,0), chrx_table,
                           by = names(chrx_table))

    # Merge the tables
    merged_table = rbind.data.frame(autosomal_table, chrx_table) %>%
      mutate(`1KG_ID` = SNPID) %>%
      mutate(NEA = Allele1) %>%
      mutate(EA = Allele2) %>%
      mutate(EAF = AF_Allele2) %>%
      mutate(P_gc = p.value) %>%
      mutate(SE_gc = SE) %>%
      select(`1KG_ID`, SNPID,
             CHR, POS, NEA, EA,AC_Allele2, EAF,
             imputationInfo, N, BETA, SE, Tstat,
             p.value, varT, varTstar, P_gc, SE_gc)
  }else{
    merged_table = autosomal_table %>%
      mutate(`1KG_ID` = SNPID) %>%
      mutate(NEA = Allele1) %>%
      mutate(EA = Allele2) %>%
      mutate(EAF = AF_Allele2) %>%
      mutate(P_gc = p.value) %>%
      mutate(SE_gc = SE) %>%
      select(`1KG_ID`, SNPID,
             CHR, POS, NEA, EA,AC_Allele2, EAF,
             imputationInfo, N, BETA, SE, Tstat,
             p.value, varT, varTstar, P_gc, SE_gc)
  }

  cat("Writing GWASSS Table...\n")
  merged_filename = gsub("XXX", "_MERGED", autosomal_pattern)
  write.table(merged_table, file = merged_filename, sep = " ",
              quote = FALSE, row.names = FALSE)
  cat(paste0("Table saved as ", merged_filename, ".\n"))

}

args = commandArgs(trailingOnly=TRUE)

if(length(args) %in% c(2, 4)){
  if(length(args) == 2){
    merge_saige(as.character(args_vec[1]),
                as.logical(args_vec[2]))
  }else{
    merge_saige(as.character(args_vec[1]),
                as.logical(args_vec[2]),
                as.character(args_vec[3]),
                as.logical(args_vec[4]))
  }
}else{
  cat("Incorrest number of arguments supplied..\n")
}

