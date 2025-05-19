mean_plot <- function(plot_mt, x_var, y_var, col_var = NULL, panel_var = NULL, plot_title = NULL, plotxlab=NULL, plotylab=NULL, scales="fixed", bar_value = 0, bar_label = c("Aveg.LSD", "Aveg.SED"), level = 0.05, basesz = 12, lineplot = TRUE){

  plot_mt <- na.omit(plot_mt)
  yMin <- min(plot_mt[, y_var])
  yMax <- max(plot_mt[, y_var])
  offSet <- 0.25 * (yMax - yMin)
  plotxlab <- ifelse(is.null(plotxlab) || any(c("NULL", "") == plotxlab), x_var, plotxlab)
  plotylab <- ifelse(is.null(plotylab) || any(c("NULL", "") == plotylab), y_var, plotylab)
  bar_label <- match.arg(bar_label)
  barMin <- yMin
  barMax <- barMin + bar_value

  if (bar_value > 0) {
    barMin <- yMin
    barMax <- barMin + bar_value
    plot_mt[, x_var] <- factor(plot_mt[, x_var], levels = c(bar_label, levels(plot_mt[, x_var])))
  }

  if (!is.null(x_var) && is.null(col_var) && is.null(panel_var)) {
    if (is.null(plot_title) || any(c("NULL", "") == plot_title)) {
      if (bar_value > 0) {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "' with ", bar_label, " (", level * 100, "%) Bar", sep="")
      } else {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "'\n", sep="")
      }
    }
    meanPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group=1))+
      labs(title=paste(plot_title, "\n", sep=""), x=paste("\n", plotxlab, sep=""), y=paste(plotylab, "\n", sep=""))+
      lims(x= levels(plot_mt[, x_var]), y = c(yMin - offSet, max(yMax + offSet, yMin + bar_value + offSet))) +
      geom_point(colour="red", size=2)+
      geom_errorbar(aes(ymax=barMax, ymin=barMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue") +
      theme_bw(basesz)
    if (lineplot) {
      meanPlot <- meanPlot + geom_line(linewidth=0.5)
    }
  }

  if (!is.null(x_var) && !is.null(col_var) && is.null(panel_var)) {
    if (is.null(plot_title) || any(c("NULL", "") == plot_title)) {
      if (bar_value > 0) {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "' by '", col_var, "' with ", bar_label, " (", level * 100, "%) Bar", sep = "")
      } else {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "' by '", col_var, "'", sep = "")
      }
    }

    meanPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]], group=.data[[col_var]], col=.data[[col_var]]))+
      labs(title=paste(plot_title, "\n", sep=""), x=paste("\n", plotxlab, sep=""), y=paste(plotylab, "\n", sep=""))+
      lims(x= levels(plot_mt[, x_var]), y = c(yMin - offSet, max(yMax + offSet, yMin + bar_value + offSet))) +
      geom_point(size=2)+
      geom_errorbar(aes(ymax=barMax, ymin=barMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
      guides(col = guide_legend(title = col_var))+
      theme_bw(basesz)
    if (lineplot) {
      meanPlot <- meanPlot + geom_line(aes(linetype=.data[[col_var]],
                                           col=.data[[col_var]]), linewidth=0.96)+
        guides(linetype = guide_legend(title = col_var))
    }
  }

  if (!is.null(x_var) && !is.null(col_var) && !is.null(panel_var)) {
    if (is.null(plot_title) || any(c("NULL", "") == plot_title)) {
      if (bar_value > 0) {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "' by '", col_var, "' for each '",
                            panel_var, "'\n with ", bar_label, " (",
                            level * 100, "%) Bar\n", sep = "")
      } else {
        plot_title <- paste("Predicted means for '", plotylab, "' vs '", plotxlab, "' by '", col_var, "' for each '", panel_var, "'\n",  sep = "")
      }
    }

    meanPlot <- ggplot(plot_mt, aes(.data[[x_var]], .data[[y_var]],
                                    group=factor(.data[[col_var]]),
                                    col=factor(.data[[col_var]])))+
      labs(title=paste(plot_title, "\n", sep=""), x=paste("\n", plotxlab, sep=""), y=paste(plotylab, "\n", sep=""))+
      lims(x= levels(plot_mt[, x_var]), y = c(yMin-0.5*offSet, max(yMax+0.5*offSet, yMin+bar_value+0.5*offSet))) +
      geom_errorbar(aes(ymax=barMax, ymin=barMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
      geom_point(size=2)+
      facet_grid(as.formula(paste("~", panel_var)), scales = scales)+
      guides(group = guide_legend(col_var))+
      guides(col = guide_legend(col_var))+
      theme_bw(basesz)
    if (lineplot) {
      meanPlot <- meanPlot + geom_line(aes(linetype=.data[[col_var]],
                                           col=.data[[col_var]]), linewidth=0.8) +
        guides(linetype = guide_legend(title = col_var))
    }

  }
  return(meanPlot)
}


