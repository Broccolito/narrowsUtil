
library(narrowsUtil)
stats = get_senlinplot_stats(par_path = par_file,
                             outfile_path = meta_file,
                             filename = stats_filename,
                             gwas_snp_name = snp_name,
                             gwas_identifier_column_name = "1KG_ID",
                             gwas_effectsize_column_name = "BETA",
                             gwas_stderr_column_name = "SE",
                             gwas_pvalue_column_name = "p.value",
                             gwas_samplesize_column_name = "N",
                             meta_snp_name = snp_name,
                             meta_identifier_column_name = "MarkerName",
                             meta_effectsize_column_name = "Effect",
                             meta_stderr_column_name = "StdErr",
                             meta_pvalue_column_name = "P-value",
                             meta_samplesize_column_name = "N",
                             check_inverse_ref = TRUE)
make_senlinplot_backend(senlinplot_stats_path = stats_filename,
                        senlinplot_filename = senlinplot_filename,
                        use_user_theme = use_user_theme,
                        user_theme = user_theme)

