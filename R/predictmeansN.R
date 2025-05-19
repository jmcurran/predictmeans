predictmeansN <- function (model, modelterm, data=NULL, pairwise=FALSE, atvar=NULL, adj="none", Df=NULL, lsd_bar=TRUE, level=0.05, covariate=NULL, meandecr=NULL, letterCI=FALSE, trans = NULL, transOff=0,responsen=NULL, count=FALSE, plotord=NULL, lineplot=TRUE, mplot=TRUE, barplot=FALSE, pplot=TRUE, plot=TRUE, prtnum=TRUE, prtplt=TRUE, permlist=NULL, ncore=3L, ndecimal=4L) {

  options(scipen=6L)
  if (any(missing(model), missing(modelterm))) {
    stop("The arguments 'model', and 'modelterm' must be provided!")
  }
  if (!(modelterm  %in%  attr(terms(model), "term.labels"))) {
    stop(paste("The", modelterm, "must be exactly a term in the model (especially check the order of interaction)."))
  }
  vars <- unlist(strsplit(modelterm, "\\:"))
  model_frame <- model.frame(model)
  if (any(!base::is.element(vars, names(model_frame)[sapply(model_frame, is.factor)]))) {
    stop(paste(vars, "must be factor(s)!"))
  }

  if (inherits(model, "aovlist")) {
    cat("Refit the 'aov' model using function 'lmer' ... \n")
    model <- aovlist_lmer(model)
  }

  if (class(model)[1] %in% c("glm", "negbin", "glmerMod", "glmmTMB")) {
    trans <- switch(class(model)[1],
                    "glm"=model$family$linkinv,
                    "negbin" = model$family$linkinv,
                    "glmerMod"=slot(model, "resp")$family$linkinv,
                    "glmmTMB"=model$modelInfo$family$linkinv)

    if (identical(trans, make.link("identity")$linkinv)) {
      trans <- NULL
    }
  }

  if (length(vars)==1) {
    atvar <- NULL
  }
  if (!is.null(permlist)) {
    pairwise <- TRUE
    if (adj=="tukey") {
      stop("The p-value can't be adjusted by Tukey methd!")
    }
  }
  if (!is.null(atvar)) {
    pairwise <- TRUE
    resvar <- vars[!(vars %in% atvar)]   # The rest of vars rather than at var
    vars <- c(atvar, resvar)
  }
  if (adj != "none") {
    pairwise <- TRUE
  }
  if (letterCI) {
    pairwise <- TRUE
    ci_level <- round(2*(1-pnorm(sqrt(2)*qnorm(1-level/2)/2)), 3)
    pmlevel <- level # <- 0.05  # back to 0.05 for comparison?
  } else {
    pmlevel <- ci_level <- level
  }

  if (!is.null(plotord)) {
    plot <- mplot <- TRUE
    if (length(plotord) != length(vars)) {
      stop(paste("plotord must be a vector of length", length(vars)))
    }
    if (is.character(plotord)) {
      i = match(plotord, vars)
      if (any(is.na(i))) {
        stop(paste("The entries", paste0(plotord[is.na(i)], sep = ", "), "do not match any term in ", modelterm))
      } else {
        plotord = i
      }
    }
  }

  if (!is.logical(meandecr)) {
    meandecr <- NULL
  } else {
    atvar <- NULL
  }

  ctr.matrix <- Kmatrix(model, modelterm, covariate, data=data, prtnum=prtnum)
  KK <- ctr.matrix$K
  label <- ctr.matrix$fctnames
  rownames(label) <- rownames(KK)
  response <- ctr.matrix$response
  mp <- mymodelparm(model)
  n.table <- table(model_frame[, vars, drop = FALSE])
  ndf <- data.frame(n.table)       ## To obtain info from model

  if (inherits(model, "glmmTMB")) {
    K <- KK
  } else {
    K <- KK[, mp$estimable, drop = FALSE]          # To match coef names
  }

  if (any(ndf$Freq == 0)) {
    rnTrt <- do.call("paste", c(ndf[, vars, drop=FALSE], sep=":"))
    rnTrt <- rnTrt[ndf$Freq!=0]      # To delete any missing level in factor
    K <- K[rownames(K) %in% rnTrt,]
    label <- label[rownames(label) %in% rnTrt,]
  }

  pm <- K %*% mp$coef
  vcovm <- mp$vcov
  ses <- as.numeric(apply(K, 1, function(x) {y <- matrix(x, nrow=1);sqrt(y %*% tcrossprod(mp$vcov, y))}))
  meanTable <- data.frame(pm, ses, label)

  if (length(Df) == 0) {
    if (inherits(model, "lme")) {
      Df <- terms(model$fixDF)[modelterm]
      names(Df) <- NULL
      meanTable$Df <- round(Df, 2)
      pairDf <- Df
      Df_diff <- Df
    } else if (inherits(model, "lmerMod")) {
      L_term <- get_contrasts_type1(model)[[modelterm]]
      Df <- mean(df_term(model, ctrmatrix=L_term))
      if (is.nan(Df)) {
        stop("You need provide Df for this model term!")
      }
      terms_df <- round(df_term(model, ctrmatrix=K), 2)
      terms_df <- terms_df[!is.na(terms_df)]
      terms_df[terms_df <= 0] <- Df
      meanTable$Df <- terms_df
    } else {
      Df <- mp$df
      meanTable$Df <- round(Df, 2)
      pairDf <- Df
      Df_diff <- Df
    }

    if (Df == 0) {
      stop("You need provide Df for this model term!")
    }
  } else {
    meanTable$Df <- Df
    pairDf <- Df
    Df_diff <- Df
  }

  if (is.null(permlist)) {
    meanTable$LL <- meanTable$pm - qt(1 - ci_level/2, df = meanTable$Df) * meanTable$ses
    meanTable$UL <- meanTable$pm + qt(1 - ci_level/2, df = meanTable$Df) * meanTable$ses
  } else {
    meanTable$LL <- meanTable$pm - 2 * meanTable$ses
    meanTable$UL <- meanTable$pm + 2 * meanTable$ses
    pmlevel <- ci_level <- level <- 0.05
    meanTable$Df <- NA
  }

  rownames(meanTable) <- NULL
  meanTable <- meanTable[c(vars, "pm", "ses", "Df", "LL", "UL")]
  meanTable[c("pm", "LL", "UL")] <- round(meanTable[c("pm", "LL", "UL")], ndecimal)
  meanTable$ses <- round(meanTable$ses, ndecimal+1)
  colnames(meanTable) <- c(vars, "Mean", "SE", "Df", paste("LL(", (1 - ci_level) * 100, "%)", sep = ""),
                           paste("UL(", (1 - ci_level) * 100, "%)", sep = ""))
  meanTable <- meanTable[do.call(order, meanTable[vars]), ]

  if (!is.null(meandecr) && is.logical(meandecr)) {
    meanTable <- meanTable[order(meanTable[["Mean"]], decreasing = meandecr), ]
  }

  if (!(is.null(trans) || identical(trans, I))) {
    meanTable[["Bk_Mean"]] <- trans(meanTable[["Mean"]]) - transOff
    meanTable[[paste("Bk_LL(", (1 - ci_level) * 100, "%)", sep = "")]] <- trans(meanTable[[paste("LL(", (1 - ci_level) * 100, "%)", sep = "")]]) - transOff
    meanTable[[paste("Bk_UL(", (1 - ci_level) * 100, "%)", sep = "")]] <- trans(meanTable[[paste("UL(", (1 - ci_level) * 100, "%)", sep = "")]]) - transOff

    if (any(names(model_frame) == response)) {    ## Transformed y before modelling
      if (inherits(model_frame[, response], "factor")) {
        model_frame$Bk_response <- as.numeric(model_frame[, response])-1
      } else {
        if (identical(trans, I) || (inherits(model, "glm") || inherits(model, "glmerMod") || inherits(model, "glmmTMB"))) {
          model_frame$Bk_response <- model_frame[, response]
          if (!is.null(dim(model_frame[, response]))) {
            model_frame$Bk_response <- model_frame[, response][,1]/rowSums(model_frame[, response])
          }
        } else {
          model_frame$Bk_response <- trans(model_frame[, response])-transOff
        }# end of if glm or glmer
      }# end of if factor
    } else {       ## Transformed y within modelling
      nresponse <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
      if (!(any(names(model_frame) == nresponse))) {
        if (is.null(responsen)) {
          stop(paste("Please provide suitable name for response variable using option responsen='", names(model_frame)[1], "'!", sep=""))
        }
        nresponse <- responsen
      }
      model_frame$Bk_response <- model_frame[, nresponse]
    }

    if (is.list(model_frame$Bk_response)) {
      model_frame$Bk_response <- model_frame$Bk_response[[1]]
    }
  } else {
    model_frame$Bk_response <- model_frame[, response]
  }

  if (letterCI) {
    if (is.null(atvar)) {
      meanTable$LetterGrp <- ci_mcp(meanTable[, grepl("^LL", names(meanTable))], meanTable[, grepl("^UL", names(meanTable))])
      meanTable <- list(Table=meanTable,
                        Note=paste("Letter-based representation of pairwise comparisons based on (",
                                   (1 - ci_level) * 100, "%) CI at ", level, " significant level.", sep=""))
    } else {
      meanTable$LetterGrp <- unlist(lapply(split(meanTable, meanTable[, atvar]),
                                           function (x) {
                                             ci_mcp(x[,grepl("^LL", names(meanTable))], x[,grepl("^UL", names(meanTable))])
                                           }))
      meanTable <- list(Table=meanTable,
                        Note=paste("Letter-based representation of pairwise comparisons based on (",
                                   (1 - ci_level) * 100, "%) CI at ", level, " significant level for each ", atvar, " group.", sep=""))
    }
  }

  if (plot) {
    if (length(vars) > 3) {
      cat("\n", "There is no plot for more than three-way interaction! \n\n")
    }

    if (is.data.frame(meanTable)) {
      plotmt <- meanTable
    } else {
      plotmt <- meanTable[[1]]
    }

    ciPlot <- ci_plot(plotmt, model_frame, plotxlab="Bk_response")
    if (prtplt) {
      print(ciPlot)
    }

    if (length(vars) == 1) {
      if (mplot) {
        meanPlot <- mean_plot(plotmt, vars, y_var="Mean", plotylab=response, bar_value = 0, level = level, lineplot = lineplot)
        if (prtplt) {
          print(meanPlot)
        }
      } # end of if mplot
      if (barplot) {
        barPlot <- bar_plot(plotmt, vars, y_var="Mean", se="SE", col_var = NULL, plotylab=response)
        if (prtplt) {
          print(barPlot)
        }
      }
    }
    if (length(vars) == 2) {
      if (is.null(plotord)) {
        plotord <- 1:2
        if (!is.null(atvar)) {
          atvar_num <- which(vars %in% atvar)
          plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
        }
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]

      if (mplot) {
        meanPlot <- mean_plot(plotmt, fact1, y_var="Mean", col_var = fact2, plotylab=response, bar_value = 0, level = level, lineplot = lineplot)
        if (prtplt) {
          print(meanPlot)
        }
      } # end if mplot
      if (barplot) {
        barPlot <- bar_plot(plotmt, fact1, y_var="Mean", se="SE", col_var = fact2, plotylab=response)
        if (prtplt) {
          print(barPlot)
        }
      }
    }

    if (all(length(vars) == 3)) {
      if (is.null(plotord)) {
        plotord <- 1:3
        if (length(atvar)==1) {
          atvar_num <- which(vars  %in%  atvar)
          plotord <- c(atvar_num, setdiff(1:length(vars), atvar_num))
        }
      }
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      fact3 <- (vars[plotord])[3]
      if (mplot) {
        meanPlot <- mean_plot(plotmt, fact1, y_var="Mean", col_var = fact2, panel_var = fact3, plotylab=response, bar_value = 0, level = level, lineplot = lineplot)
        if (prtplt) {
          print(meanPlot)
        }
      }

      if (barplot) {
        barPlot <- bar_plot(plotmt, fact1, y_var="Mean", se="SE", col_var = fact2, panel_var = fact3, plotylab=response)
        if (prtplt) {
          print(barPlot)
        }
      }
    }
  }

  outputlist <- list(mean_table=meanTable, predictmeansPlot=meanPlot, predictmeansBarPlot=barPlot,
                     predictmeansciPlot=ciPlot)
  class(outputlist) <- "pdmlist"
  return(outputlist)

}
