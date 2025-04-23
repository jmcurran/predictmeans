#' Plot the predicted means
#'
#' Plot the predicted means
#'
#'
#' @param x an object of class \code{pdmlist}. See \link{list("predictmeans")}
#' for more info.
#' @param type Plot type. Can be one of \code{"bar"}, \code{"mean"},
#' \code{"pm"}. One letter abbrevitions, i.e. \code{"b"}, \code{"m"},
#' \code{"pm"} will work as well.
#' @param ... other arguments
#' @importFrom ggplot2 ggplot aes geom_bar
#' @export plot.pdmlist
plot.pdmlist = function(x, type = c("bar", "mean", "pm"), ...){
  if (length(vars) > 3){
    cat("\n", "There is no plot for more than three-way interactions! \n\n")
  }

  mt = x$meanTable
  atvar = x$atvar

  plotmt <- na.omit(mt)
  yMin <- min(plotmt[, "pm"])
  yMax <- max(plotmt[, "pm"])
  offSet <- 0.25 * (yMax - yMin)

  if (is.null(atvar) || all(atvar%in%c("NULL", ""))) {
    if (lsd_bar){
      LSD_value <- LSD[3]
    }else{
      LSD_value <- SED.out[3]
    }
  }else{
    if (lsd_bar){
      LSD_value <- mean(attr(LSD,"For the Same Level of Factor")[1, atvar])
    }else{
      LSD_value <- mean(attr(SED.out, "For the Same Level of Factor")[1, atvar])
    }
  }

  if (lsd_bar){
    bar_label <- "Aveg.LSD"
  }else{
    bar_label <- "Aveg.SED"
  }

  up <- yMin + LSD_value
  lsdBar <- cbind(plotmt[, vars, drop=FALSE], up=up, yMin=yMin)
  limits <- aes(ymax = (pm + ses)*(pm > 0) + pmin(pm + ses, 0)*(pm <= 0), ymin=(pm - ses)*(pm < 0) + pmax(pm - ses, 0)*(pm >= 0))

  if (length(vars) == 1) {
    mxlab <- plotxlab

    if (is.null(plotxlab) || plotxlab%in%c("NULL", "")){
      mxlab <- paste("\n", vars, sep="")
    }

    mylab <- plotylab

    if (is.null(plotylab) || plotylab%in%c("NULL", "")){
      mylab <- paste(response, "\n", sep="")
    }

    if (mplot) {
      if (newwd){
        dev.new()
      }

      mtitle <- plottitle

      if (is.null(plottitle) || plottitle%in%c("NULL", "")){
        mtitle <- paste("Predicted means for \"", vars, "\" with ", bar_label, " (", slevel * 100, "%) Bar", sep="")
      }

      p1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm, group=1))+
        labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
        lims(x= c(bar_label, levels(plotmt[, vars])), y = c(yMin - offSet, max(yMax + offSet, yMin + LSD_value + offSet))) +
        geom_point(colour="red", size=2)+
        #  geom_line(linewidth=0.5)+
        geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue") +  #data=lsdBar,
        theme_bw(basesz)

      if (lineplot){
        p1 <- p1 + geom_line(linewidth=0.5)
      }

      meanPlot <- p1

      if (prtplt){
        print(p1)
      }
    } # end of if mplot

    if (barplot) {

      if (newwd){
        dev.new()
      }

      mtitle <- plottitle

      if (is.null(plottitle) || plottitle%in%c("NULL", "")){
        mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
      }

      dodge <- position_dodge(width=0.9)
      bp1 <- ggplot(plotmt, aes(eval(parse(text = vars)), pm))+
        labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
        ylim(c(min(0, (min(pm) - max(ses)) * 1.2), max(0, (max(pm) + max(ses)) * 1.2))) +
        geom_bar(position=dodge, stat="identity", fill="papayawhip", colour="darkgreen") +
        geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())

      predictmeansBarPlot <- bp1
      if (prtplt){
        print(bp1)
      }
    }
  } else if (length(vars) == 2) {
    ##NB: I have changed this to an else statement because there shouldn't be any other way to get here. That is, length(vars) will not
    ##change from 1 to 2 between the block where length(vars)==1 is true and here, therefore the cases should be else if statements.
    if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
      plotord <- 1:2
      if (!(is.null(atvar) || all(atvar%in%c("NULL", "")))) {
        atvar_num <- which(vars %in% atvar)
        plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
      }
    }

    fact1 <- (vars[plotord])[1]
    fact2 <- (vars[plotord])[2]

    mxlab <- plotxlab

    if (is.null(plotxlab) || plotxlab%in%c("NULL", "")) {
      mxlab <- paste("\n", fact1, sep="")
    }

    mylab <- plotylab

    if (is.null(plotylab) || plotylab%in%c("NULL", "")) {
      mylab <- paste(response, "\n", sep="")
    }

    if (mplot) {

      if (newwd){
        dev.new()
      }

      mtitle <- plottitle

      if (is.null(plottitle) || plottitle%in%c("NULL", "")){
        mtitle <- paste("Predicted means for \"", fact1, "\" by \"", fact2, "\" with ",
                        paste(atvar, collapse=" "), " ", bar_label, " (", slevel * 100, "%) Bar", sep = "")
      }

      plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))
      p2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), col=eval(parse(text = fact2))))+
        labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
        lims(x= levels(plotmt[, fact1]), y = c(yMin - offSet, max(yMax + offSet, yMin + LSD_value + offSet))) +
        geom_point(size=2)+
        # geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), linewidth=0.96)+
        geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
        #  guides(linetype = guide_legend(title = fact2))+
        guides(col = guide_legend(title = fact2))+
        theme_bw(basesz)

      if (lineplot) {
        p2 <- p2 + geom_line(aes(linetype=eval(parse(text = fact2)),
                                 col=eval(parse(text = fact2))), linewidth=0.96)+
          guides(linetype = guide_legend(title = fact2))
      }

      meanPlot <- p2

      if (prtplt){
        print(p2)
      }
    } # end if mplot

    if (barplot) {

      if (newwd) {
        dev.new()
      }

      mtitle <- plottitle

      if (is.null(plottitle) || plottitle%in%c("NULL", "")) {
        mtitle <- paste("Predicted means for \"", modelterm, "\" with SE Bars", sep = "")
      }

      dodge <- position_dodge(width=0.9)

      bp2 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=eval(parse(text = fact2)), fill=  eval(parse(text = fact2))))+
        geom_point(position=dodge) +
        geom_bar(stat = "identity", position=dodge)+
        labs( x = mxlab, y = mylab, title =mtitle)+
        ylim(c(min(0, (min(pm) - max(ses)) * 1.1), max(0, (max(pm) + max(ses)) * 1.1))) +
        geom_errorbar(limits, position=dodge, width=0.25, colour="blue")+theme_bw(basesz)+
        scale_fill_brewer()+
        theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank())+
        theme(legend.position = "top")+guides(fill = guide_legend(title = fact2))
      predictmeansBarPlot <- bp2

      if (prtplt) {
        print(bp2)
      }
    }
  } else if (all(length(vars) == 3, mplot)) {
    ## NB: I have changed this to else if, because neither length(vars)
    ## or mplot should have changed, and I believe this is just the third
    ## choice
    if (newwd) {
      dev.new()
    }

    mtitle <- plottitle

    if (is.null(plotord) || all(plotord%in%c("NULL", ""))) {
      plotord <- 1:3

      if (!(is.null(atvar) || all(atvar%in%c("NULL", ""))) && length(atvar)==1){
        atvar_num <- which(vars %in% atvar)
        plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
      }
    }

    fact1 <- (vars[plotord])[1]
    fact2 <- (vars[plotord])[2]
    fact3 <- (vars[plotord])[3]
    plotmt[, fact1] <- factor(plotmt[, fact1], levels = c(bar_label, levels(plotmt[, fact1])))

    if (is.null(plottitle) || plottitle%in%c("NULL", "")){
      mtitle <- paste("Predicted means for '", fact1, "' by '", fact2, "' for each '",
                      fact3, "'\n with ", paste(atvar, collapse=" "), " ",
                      bar_label, " (", slevel * 100, "%) Bar\n", sep = "")
    }

    mxlab <- plotxlab

    if (is.null(plotxlab) || plotxlab%in%c("NULL", "")){
      mxlab <- paste("\n", fact1, sep="")
    }

    mylab <- plotylab

    if (is.null(plotylab) || plotylab%in%c("NULL", "")) {
      mylab <- paste(response, "\n", sep="")
    }

    p3 <- ggplot(plotmt, aes(eval(parse(text = fact1)), pm, group=factor(eval(parse(text = fact2))), col=factor(eval(parse(text = fact2)))))+
      labs(title=paste(mtitle, "\n", sep=""), x=mxlab, y=mylab)+
      lims(x= levels(plotmt[, fact1]), y = c(yMin-0.5*offSet, max(yMax+0.5*offSet, yMin+LSD_value+0.5*offSet))) +
      geom_errorbar(aes(ymax=up, ymin=yMin, x=bar_label), width=0.15, linewidth=0.8, colour="blue")+
      geom_point(size=2)+
      # geom_line(aes(linetype=eval(parse(text = fact2)), col=eval(parse(text = fact2))), linewidth=0.8)+
      facet_grid(eval(parse(text = paste("~",fact3, sep=""))))+
      guides(group = guide_legend(fact2))+
      # guides(linetype = guide_legend(fact2))+
      guides(col = guide_legend(fact2))+
      theme_bw(basesz)

    if (lineplot){
      p3 <- p3 + geom_line(aes(linetype=eval(parse(text = fact2)),
                               col=eval(parse(text = fact2))), linewidth=0.8) +
        guides(linetype = guide_legend(title = fact2))
    }

    meanPlot <- p3

    if (prtplt){
      print(p3)
    }
  }
}
