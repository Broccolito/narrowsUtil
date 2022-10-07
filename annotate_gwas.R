if(!require("BioThingsClient.R")){
  install_github("biothings/BioThingsClient.R")
  library("BioThingsClient.R")
}

if(!require("devtools")){
  install.packages("devtools")
  library("devtools")
}

if(!require("dplyr")){
  install.packages("dplyr")
  library("dplyr")
}

if(!require("data.table")){
  install.packages("data.table")
  library("data.table")
}

if(!require("purrr")){
  install.packages("purrr")
  library("purrr")
}

if(!require("ggplot2")){
  install.packages("ggplot2")
  library("ggplot2")
}

if(!require("writexl")){
  install.packages("writexl")
  library("writexl")
}

load_data = function(file_name = "FHS_EA_MRS_5e8_snplist.txt",
                     by_marker = TRUE,
                     marker_name_column = "SNPID",
                     chr_column = "CHR",
                     pos_column = "POS",
                     ref_column = "NEA",
                     alt_column = "EA",
                     pvalue_column = "p.value"){
  
  cat("Loading GWAS Statistics...\n")
  raw_data = fread(file_name) %>%
    as.data.frame()
  
  if(by_marker){
    marker_name = as.character(raw_data[[marker_name_column]])
    marker_name = strsplit(marker_name, split = "[:]")
    chr = unlist(lapply(marker_name, function(x){x[1]}))
    pos = unlist(lapply(marker_name, function(x){x[2]}))
    ref = unlist(lapply(marker_name, function(x){x[3]}))
    alt = unlist(lapply(marker_name, function(x){x[4]}))
  }else{
    chr = raw_data[[chr_column]]
    pos = raw_data[[pos_column]]
    ref = raw_data[[ref_column]]
    alt = raw_data[[alt_column]]
  }
  
  marker = raw_data[[marker_name_column]]
  chr[chr=="X"] = 23
  chr = as.numeric(chr)
  pos = as.numeric(pos)
  pvalue = as.numeric(raw_data[[pvalue_column]])
  query = paste0("chr", chr, ":g.", pos, ref, ">", alt)
  
  if(by_marker){
    raw_data = raw_data %>%
      select(-all_of(marker_name_column),
             -all_of(pvalue_column))
  }else{
    raw_data = raw_data %>%
      select(-all_of(marker_name_column),
             -all_of(chr_column),
             -all_of(pos_column),
             -all_of(ref_column),
             -all_of(alt_column),
             -all_of(pvalue_column))
  }
  
  data = tibble(
    marker = marker,
    chr = chr,
    pos = pos,
    ref = ref,
    alt = alt,
    pvalue = pvalue,
    query = query
  ) %>%
    cbind.data.frame(raw_data)
  
  cat("GWAS Statistics Loaded...\n")
  return(data)
}

annotate_data = function(data, file_name = "data_annotated.xlsx"){
  variant_client = BioThingsClient("variant")
  annotation = vector()
  for(q in data$query){
    try({
      v = btGet(variant_client, q)[[1]]
      
      gene_symbol = NA
      gene_name = NA
      variant_type = NA
      variant_consequence = NA
      variant_CADD_phred = NA
      variant_phylop_vertebrate = NA
      variant_polyphen = NA
      variant_sift = NA
      variant_rs = NA
      
      gene_symbol = paste(v[["dbsnp"]][["gene"]][["symbol"]], collapse = "; ")
      gene_name = gsub(pattern = ",", ";", v[["dbsnp"]][["gene"]][["name"]])
      variant_consequence = gsub(",","; ", paste(v[["cadd"]][["consequence"]],
                                                 collapse = "; "))
      variant_CADD_phred = v[["cadd"]][["phred"]]
      variant_phylop_vertebrate = v[["cadd"]][["phylop"]][["vertebrate"]]
      variant_polyphen = v[["cadd"]][["polyphen"]][["val"]]
      variant_sift = v[["cadd"]][["sift"]][["val"]]
      variant_rs = v[["dbsnp"]][["rsid"]]
      
      gene_symbol = ifelse(is.null(gene_symbol), NA, gene_symbol)
      gene_name = ifelse(is.null(gene_name), NA, gene_name)
      variant_consequence = ifelse(is.null(variant_consequence), NA, variant_consequence)
      
      variant_CADD_phred = ifelse(is.null(variant_CADD_phred), NA, variant_CADD_phred)
      variant_phylop_vertebrate = ifelse(is.null(variant_phylop_vertebrate),
                                         NA, variant_phylop_vertebrate)
      variant_polyphen = ifelse(is.null(variant_polyphen), NA, variant_polyphen)
      variant_sift = ifelse(is.null(variant_sift), NA, variant_sift)
      variant_rs = ifelse(is.null(variant_rs), NA, variant_rs)
      
      anno = tibble(query = q, variant_rs, gene_symbol, gene_name,
                    variant_consequence, variant_CADD_phred, 
                    variant_phylop_vertebrate, variant_polyphen,
                    variant_sift)
      
      annotation = rbind.data.frame(annotation, anno)
      cat(paste0("Annotating ", q," ..\n"))
    }, silent = TRUE)
  }
  
  data_annotated = data %>%
    left_join(annotation, by = "query")
  
  write_xlsx(data_annotated, path = file_name)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) %in% c(2, 4)){
  if(length(args) == 2){
    data = load_data(file_name = as.character(args[1]))
    annotate_data(data, args[2])
  }else{
    data = load_data(file_name = as.character(args[1]), 
                     marker_name_column = as.character(args[3]),
                     pvalue_column = as.character(args[4]))
    annotate_data(data, args[2])
  }
}else{
  cat("Incorrest number of arguments supplied..\n")
}
