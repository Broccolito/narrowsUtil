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

if(!require("latex2exp")){
  install.packages("latex2exp")
  library("latex2exp")
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

plot_manhattan = function(data,
                          color1 = "cornflowerblue",
                          color2 = "blue4",
                          threshold_linecolor = "black",
                          plot_name = "manhattan_plot.png"){

  cat("Reshaping data to Manhattan plot data frame...\n")
  data_chromosome_position = data %>%
    group_by(chr) %>%
    summarise(chromosome_length = max(pos)) %>%
    mutate(tot = cumsum(as.numeric(chromosome_length)) - chromosome_length) %>%
    select(-chromosome_length)

  data_cleaned = left_join(data, data_chromosome_position, by = "chr") %>%
    arrange(chr, pos) %>%
    mutate(bp_cum = pos + tot) %>%
    mutate(nlogp = -log(pvalue, base = 10)) %>%
    select(chr, pos, pvalue, nlogp, tot, bp_cum)

  axisdf = data_cleaned %>%
    group_by(chr) %>% summarize(center = (max(bp_cum) + min(bp_cum))/2)

  plt = ggplot(data_cleaned, aes(x = bp_cum, y = nlogp)) +
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c(color1, color2), 22)) +
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(limits = c(0, max(max(data_cleaned$nlogp)+1,
                                         -log(5e-8,base = 10)+1)),
                       expand = c(0, 0)) +
    geom_hline(yintercept = -log(5e-8,base = 10),
               linetype = "solid", color = threshold_linecolor) +
    geom_hline(yintercept = -log(1e-5,base = 10),
               linetype = "dashed", color = threshold_linecolor) +
    xlab("Chromosome") +
    ylab(latex2exp::TeX("-$log_{10}$(p-value)")) +
    theme_pubr() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20)
    )

  cat("Saving the Manhattan Plot as a PNG file...\n")
  ggsave(filename = plot_name, plot = plt, device = "png", dpi = 1200,
         height = 6, width = 15)
  cat(paste0("Plot saved as ", plot_name, "\n"))
}

plot_qq = function(data, plot_name = "qq_plot.png"){

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
    xlab(latex2exp::TeX("Expected $-log_{10}$(p-value)")) +
    ylab(latex2exp::TeX("Observed $-log_{10}$(p-value)")) +
    ggtitle(label = NULL, subtitle = TeX(paste0("$\\lambda = $", lambda))) +
    theme_classic() +
    theme(text = element_text(size = 15))

  cat("Saving the Q-Q plot as png...\n")
  ggsave(filename = plot_name, device = "png", plot = plt,
         dpi = 1200, width = 6, height = 4.5)
  cat(paste0("Plot saved as ", plot_name, "\n"))
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) %in% c(2, 4)){
  if(length(args) == 2){
    merged_file = as.character(args[1])
    plot_name_suffix = as.character(args[2])
    manhattan_plot_name = paste0(plot_name_suffix, "_manhattan.png")
    qq_plot_name = paste0(plot_name_suffix, "_qq.png")
    cat("Loading GWAS summary statistics...\n")
    data = load_data(file_name = merged_file)
    cat("Making the Manhattan plot...\n")
    plot_manhattan(data = data, plot_name = manhattan_plot_name)
    cat("Make the Q-Q plot...\n")
    plot_qq(data = data, plot_name = qq_plot_name)
  }else{
    merged_file = as.character(args[1])
    plot_name_suffix = as.character(args[2])
    markername = as.character(args[3])
    pvalue = as.character(args[4])
    manhattan_plot_name = paste0(plot_name_suffix, "_manhattan.png")
    qq_plot_name = paste0(plot_name_suffix, "_qq.png")
    cat("Loading GWAS summary statistics...\n")
    data = load_data(file_name = merged_file,
                     marker_name_column = markername,
                     pvalue_column = pvalue)
    cat("Making the Manhattan plot...\n")
    plot_manhattan(data = data, plot_name = manhattan_plot_name)
    cat("Make the Q-Q plot...\n")
    plot_qq(data = data, plot_name = qq_plot_name)
  }
}else{
  cat("Incorrest number of arguments supplied..\n")
}


