#' Make a Q-Q Plot
#'
#' This function makes a Q-Q plot with the default data frame of GWAS
#' summary statistics from load_data function
#'
#' @param data an data frame object loaded from load_data function
#' @param plot_name name of the Q-Q plot file
#' @return NULL
#' @export
make_qqplot_backend = function(data,
                               plot_name = "qq_plot.png"){

  cat("Calculating Lambda...\n")
  p_values = data[["pvalue"]]
  lambda = median(qchisq(p_values, df = 1, lower.tail = FALSE))/qchisq(0.5,df=1)
  lambda = round(lambda,2)

  cat("Sorting P Values...\n")
  n = length(p_values)
  j = 1:n
  x = (j-0.5)/n
  d = sort(p_values)

  cat("Calculating plot x,y limit... \n")
  xylim = max(c(8,round(-log10(min(p_values))+2,0)))

  cat("Calculating confidence intervals... \n")
  u95 = -log10(qbeta(0.025, j, (n-j+1)))
  l95 = -log10(qbeta(0.975, j, (n-j+1)))

  qqdata = tibble(x = -log10(x),
                  y = -log10(d),
                  u95, l95)

  cat("Generating Q-Q plot...\n")
  plt = ggplot(data = qqdata, aes(x = x, y = y)) +
    geom_point(shape = 21, color = "black", fill = "gray", size = 2) +
    geom_line(aes(y = x)) +
    geom_line(aes(y = u95), linetype = "dotted") +
    geom_line(aes(y = l95), linetype = "dotted") +
    xlim(c(0, xylim)) +
    ylim(c(0, xylim)) +
    xlab(TeX("Expected $-log_{10}$(p-value)")) +
    ylab(TeX("Observed $-log_{10}$(p-value)")) +
    ggtitle(label = NULL, subtitle = TeX(paste0("$\\lambda = $", lambda))) +
    theme_classic() +
    theme(text = element_text(size = 15))

  cat("Saving the Q-Q plot as png...\n")
  ggsave(filename = plot_name, device = "png", plot = plt,
         dpi = 1200, width = 6, height = 4.5)
  cat(paste0("Plot saved as ", plot_name, "\n"))
}
