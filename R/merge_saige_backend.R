#' @export
merge_saige_backend = function(autosomal_pattern = "FHS_EA_MRS_chrXXX.txt",
                               autosomal_only = FALSE,
                               chrx_filename = "FHS_EA_MRS_plink2_run_chrX.txt",
                               write_file = TRUE){

  file_list = sapply(1:22, function(x){
    gsub(pattern = "XXX", replacement = x, autosomal_pattern)
  }) %>%
    as.list()

  file_lookup = map(file_list, function(x){
    tibble(file = x, found = file.exists(x))
  }) %>%
    reduce(rbind.data.frame)

  print(kable(file_lookup))

  if(sum(!file_lookup$found)>0){
    stop("Some of the files were not found. Please check input...")
  }

  cat("Reading GWAS Summary Statistics...\n")
  merged_table = map(file_list, fread) %>%
    reduce(rbind.data.frame) %>%
    as.data.frame() %>%
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

  if(!autosomal_only){

    filter_chrx_table = function(chrx_filename, sample_size){
      d = read.delim(file = chrx_filename, sep = "\t")
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

    # Filter sex chromosome table
    sample_size = max(merged_table$N)
    chrx_table = filter_chrx_table(chrx_filename, sample_size)
    chrx_table = full_join(head(merged_table,0), chrx_table,
                           by = names(chrx_table))

    # Merge the tables
    merged_table = rbind.data.frame(merged_table, chrx_table) %>%
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

  cat("Saving to merged table...\n")
  merged_filename = gsub("XXX", "_MERGED", autosomal_pattern)
  write.table(merged_table, file = merged_filename, sep = " ",
              quote = FALSE, row.names = FALSE)
}
