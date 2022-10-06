#' Load GWAS Summary Statistics from file
#'
#' This function loads GWAS summary statistics from file
#'
#' @param file_name The name of the file containing the statistics
#' @param by_marker A Boolean variable indicating whether the marker names will be used to specify the variant. If TRUE, marker_name_column has to be non-empty. If FALSE, chr_column, pos_column, ref_column, alt_column have to be non-empty
#' @param marker_name_column The column name of the 1KG markers, usually in the format of chr:pos:ref:alt
#' @param chr_column, The column name of the chromosome number column
#' @param chr_column, The column name of the chromosome number column
#' @param pos_column, The column name of the chromosome position column
#' @param ref_column, The column name of the reference allele column
#' @param alt_column, The column name of the alternative allele column
#' @param pvalue_column, The column name of the p-value column
#' @param delimiter, The delimiter used in the file
#' @return data, a data frame of the summary statistics
#' @export
load_gwasss_backend = function(file_name = "gwas_summary_stats.txt",
                               by_marker = TRUE,
                               marker_name_column = "SNPID",
                               chr_column = "CHR",
                               pos_column = "POS",
                               ref_column = "NEA",
                               alt_column = "EA",
                               pvalue_column = "p.value",
                               delimiter = " "){

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

