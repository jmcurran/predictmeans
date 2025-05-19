ci_plot <- function(plot_mt, point_df, plotxlab, jitterv=0.2, basesz = 12){

  plot_mt <- na.omit(plot_mt)
  LL <- Mean_v <- Treat <- UL <- NULL
  vars <- names(plot_mt)[sapply(plot_mt, is.factor)]
  plot_mt$Treat <- do.call("paste", c(plot_mt[, vars, drop=FALSE], sep=":"))
  point_df$Treat <- do.call("paste", c(point_df[, vars, drop=FALSE], sep=":"))

  slevel <- gsub("UL", "", names(plot_mt)[grepl("^UL", names(plot_mt))])
  if (any(grepl("^Bk_", names(plot_mt)))) {
    plot_mt[["Mean_v"]] <- plot_mt[["Bk_Mean"]]
    names(plot_mt)[grepl("^Bk_UL", names(plot_mt))] <- "UL"
    names(plot_mt)[grepl("^Bk_LL", names(plot_mt))] <- "LL"
    mtitle <- paste("Plot Back Transformed Means with ", slevel, " CIs for '", paste(vars, collapse=":"), "'", "\n", sep = "")
  } else {
    plot_mt[["Mean_v"]] <- plot_mt[["Mean"]]
  names(plot_mt)[grepl("^UL", names(plot_mt))] <- "UL"
  names(plot_mt)[grepl("^LL", names(plot_mt))] <- "LL"
  mtitle <- paste("Plot Means with ", slevel, " CIs for '", paste(vars, collapse=":"), "'", "\n", sep = "")
  }

  xMax <- max(max(plot_mt[["UL"]]), point_df[[plotxlab]], na.rm=TRUE)
  xMin <- min(min(plot_mt[["LL"]]), point_df[[plotxlab]], na.rm=TRUE)
  xoffSet <- 0.15 * (xMax - xMin)
  xlimv <- c(xMin - xoffSet, xMax + xoffSet)

  ciPlot <- ggplot(plot_mt, aes(Mean_v, Treat)) +
    labs(title=mtitle, x="", y="") + # paste0("\n", plotxlab)
    xlim(xlimv) +
    geom_point(colour="red") + geom_errorbarh(aes(xmax = UL, xmin=LL), height=0.2, linewidth=0.8, colour="red") +
    scale_y_discrete(limits = rev(plot_mt[["Treat"]]))+
    geom_point(aes(x=.data[[plotxlab]], y=Treat), shape=1,
               position = position_jitter(width = jitterv,
                                          height = jitterv),
               colour="blue", alpha=0.6, data=point_df)+
    theme_bw(basesz)

  return(ciPlot)
}



