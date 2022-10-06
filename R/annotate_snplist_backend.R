#' Annotate Summary Statistics
#'
#' This function annotates GWAS summary statistics using the Biothings API
#'
#' @param data an data frame object loaded from load_data function
#' @param file_name The name of the file where the annotated results are saved
#' @param write_file A boolean value indicating whether the annotated results need to be saved
#' @return data_annotated, a data frame of the annotated summary statistics
#' @export
annotate_snplist_backend = function(data,
                                    file_name = "data_annotated.xlsx",
                                    write_file = TRUE
){
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

  return(data_annotated)
}
