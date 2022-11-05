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

make_crosstrait_senlin_plot_from_stats = function(senlinplot_stats_path = "forestplot_stats.csv",
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
               footnote = paste0("\n    ", d$marker[1]),
               theme = tm); plt
  
  p_wh = get_wh(plot = plt, unit = "in")
  ggsave(filename = senlinplot_filename, plot = plt,
         device = "png",
         dpi = 1200,
         width = p_wh[1], height = p_wh[2], units = "in")
  
}


args = commandArgs(trailingOnly=TRUE)

if(length(args) == 2){
  
  stats_filename = as.character(args[1])
  plot_name = as.character(args[2])
  
  make_crosstrait_senlin_plot_from_stats(senlinplot_stats_path = stats_filename,
                              senlinplot_filename = plot_name)
  
  if(file.exists("Rplots.pdf")){
    file.remove("Rplots.pdf")
  }
  
}else{
  cat("Incorrect number of arguments supplied..\n")
}


