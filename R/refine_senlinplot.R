refine_senlinplot = function(senlinplot_stats_path = "forestplot_stats.csv",
                             senlinplot_filename = "senlinplot.png",
                             use_user_theme = FALSE,
                             user_theme = NA){

  d = read.csv(senlinplot_stats_path)
  d$study_name = strsplit(d$path, split = "[/]") %>%
    lapply(function(x){
      x[length(x)]
    }) %>%
    unlist()
  d = d %>%
    select(path, study_name, marker, effect, se, pvalue,
           c95_lower, c95_upper, n) %>%
    filter(!is.na(effect)) %>%
    mutate(effect = signif(effect, 3)) %>%
    mutate(se = round(se, 2)) %>%
    mutate(c95_lower = round(c95_lower, 2)) %>%
    mutate(c95_upper = round(c95_upper, 2)) %>%
    mutate(ci = paste0(effect," (", c95_lower, ", ", c95_upper, ")"))

  names(d) = c("path", "Study", "marker", "Effect", "se", "P-value",
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

  if(use_user_theme){

    if(is.na(user_theme)){
      cat("No user theme supplied - use default theme...\n")
    }else{
      tm = user_theme
    }

  }

  plt = forest(data = select(d, Study, Effect, `P-value`, `95% CI`, N, ` `),
               est = d$Effect,
               lower = d$c95_lower,
               upper = d$c95_upper,
               ci_column = 6,
               is_summary = c(rep(FALSE, dim(d)[1]-1), TRUE),
               footnote = paste0("\n", d$marker[1]),
               theme = tm); plt

  p_wh = get_wh(plot = plt, unit = "in")
  ggsave(filename = senlinplot_filename, plot = plt,
         device = "png",
         dpi = 1200,
         width = p_wh[1], height = p_wh[2], units = "in")

}
