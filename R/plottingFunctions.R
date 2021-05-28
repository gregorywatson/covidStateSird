plotCols <- c("#AEC441", "#008875", "#0094D6", "#F2E000", "#C362A6", "#F78F1E", "#2E62CC", "#FFC024")

#' @export
plotCumulativeCases <- function(allS, state, days, statePop, plotT, endPlot, plotCol = plotCols[2]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qC <- t(apply(statePop - allS, 1, quantile, probs = c(.05, .5,  .95),  na.rm =     TRUE))[1:length(plotT),]
  
  plot(plotT, qC[,2], col = NA, ylim = c(0, max(qC[,3])), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", xlab = "")
  
  polygon(c(plotT, rev(plotT)),
          c(qC[,1], rev(qC[,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT, qC[,2], col = plotCol, lwd = 2)
  points(days, state[1,], pch = 19, col = "grey33")
  
  # x-axis tick marks
  if(max(qC) > 4000000) {
    axisTicks <- seq(0,10000000,1000000)
  }
  else if(max(qC) > 2000000) {
    axisTicks <- seq(0,5000000,500000)
  }
  else if(max(qC) > 1000000) {
    axisTicks <- seq(0,2000000,250000)
  } else if(max(qC) > 500000) {
    axisTicks <- seq(0,2000000,100000)
  } else if(max(qC) > 100000) {
    axisTicks <- seq(0,500000,50000)
  } else if (max(qC) > 10000) {
    axisTicks <- seq(0,500000,5000)
  } else if (max(qC) > 1000) {
    axisTicks <- seq(0,500000,500)
  } else if (max(qC) > 100) {
    axisTicks <- seq(0,500000,50)
  } else {
    axisTicks <- seq(0,500000,5)
  }
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33", col   = "grey33", cex.axis = 1.25)
  
  mtext("Cumulative Confirmed Cases", side=2, line=3.1, col="grey33", cex=1)
}

#' @export
plotActiveCases <- function(allI, state, plotT, endPlot, plotCol = plotCols[5]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qI <- t(apply(allI, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
  plot(plotT, qI[,2], col = NA, ylim = c(0, max(qI[,3]) * 1.04), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", yaxs = "i", xlab = "")
  
  polygon(c(plotT, rev(plotT)),
          c(qI[,1], rev(qI[,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT, qI[,2], col = plotCol, lwd = 2)
  
  # x-axis tick marks
  if(max(qI) > 4000000) {
    axisTicks <- seq(0,10000000,1000000)
  }
  else if(max(qI) > 2000000) {
    axisTicks <- seq(0,5000000,500000)
  }
  else if(max(qI) > 1000000) {
    axisTicks <- seq(0,2000000,250000)
  } else if(max(qI) > 500000) {
    axisTicks <- seq(0,2000000,100000)
  } else if(max(qI) > 100000) {
    axisTicks <- seq(0,500000,50000)
  } else if (max(qI) > 10000) {
    axisTicks <- seq(0,500000,5000)
  } else if (max(qI) > 1000) {
    axisTicks <- seq(0,500000,500)
  } else if (max(qI) > 100) {
    axisTicks <- seq(0,500000,50)
  } else {
    axisTicks <- seq(0,500000,5)
  }
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
  mtext("Actively Infected Confirmed Cases", side=2, line=3.1, col="grey33", cex=1)
}

#' @export
plotDailyDeaths <- function(allD, allStateFit, stateFit, state, days, plotT, endPlot, SAMPS, rfError, plotCol = plotCols[6]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  allDD <- apply(allD, 2, diff)
  qDD <- t(apply(allDD, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
  if(rfError) {
    wideDeathLo <- matrix(allStateFit$deathLo, nrow = length(stateFit$times), ncol = ncol(allD), byrow = F)
    qDD[,1] <- t(apply(apply(wideDeathLo, 2, diff), 1, quantile, probs = c(.05,.5),  na.rm = TRUE))[1:length(plotT),1]
    wideDeathUp <- matrix(allStateFit$deathUp, nrow = length(stateFit$times), ncol = ncol(allD), byrow = F)
    qDD[,3] <- t(apply(apply(wideDeathUp, 2, diff), 1, quantile, probs = c(.95,.5),  na.rm = TRUE))[1:length(plotT),1]
  }
  
  plot(plotT, qDD[,2], col = NA, ylim = c(0, max(qDD[,3]) * 1.04), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", yaxs = "i", xlab = "")
  
  polygon(c(plotT, rev(plotT)),
          c(qDD[,1], rev(qDD[,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT, qDD[,2], col = plotCol, lwd = 2)
  points(days[1:(length(days)-1)], diff(state[2,]), pch = 19, col = "grey33")
  
  if(max(qDD) > 1000) {
    axisTicks <- seq(0,100000,100)
  } else if (max(qDD) > 500) {
    axisTicks <- seq(0,100000,100)
  } else if (max(qDD) > 100) {
    axisTicks <- seq(0,100000,50)
  } else if (max(qDD) > 50) {
    axisTicks <- seq(0,100000,10)
  } else if (max(qDD) > 10) {
    axisTicks <- seq(0,100000,5)
  } else {
    axisTicks <- seq(0,100000,1)
  }
  
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
  mtext("Daily Deaths", side=2, line=3.1, col="grey33", cex=1)
}

#' @export
plotRt <- function(allRt, state, days, plotT, endPlot, plotCol = plotCols[1]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qRt <- t(apply(allRt, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
  plot(plotT, qRt[,2], col = NA, ylim = c(0, max(qRt[,3], na.rm = T)), xaxs = "i",
       xaxt = "n", yaxt = "n", ylab = "", xlab = "")
  abline(h = 1, col = "grey66", lty = 2)
  polygon(c(plotT, rev(plotT)),
          c(qRt[,1], rev(qRt[,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT, qRt[,2], col = plotCol, lwd = 2)
  points(days, state[1,], pch = 19, col = "grey33")
  
  axisTicks <- seq(0,4,.5)
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = axisTicks, cex = .9, tick = F, hadj = .75)
  
  mtext("Estimated Rt 10-day Moving Average", side=2, line=3.1, col="grey33", cex=1)
  
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
  #  mtext(stateAbbrev, outer=TRUE,  cex=2, line=-4)
}

#' @export
plotVelocityFit <- function(posteriorMean, statesLong, states, stateInterventions, fileName = NULL) {
  if(!is.null(fileName)) {
    pdf(file = fileName, width = 6, height = 6)
  }
  for(i in 1:length(states)) {
    testV <- velocitiesState(statesLong, states[i], stateInterventions, minCases = 100)
    tt <- testV$cases$t
    y <- testV$cases$y
    y[which(y <= 0)] <- NA
    yy <-(y)
    
    plot(tt,(yy), col = "white", pch = 19, type = "o", lwd = .5, cex = .7,
         main = states[i], ylim = range((yy), na.rm = T), ylab = "d/dt Log Cumulative Cases",   xlab = "Days Since 100+ Cases")
    intv <- (testV$cases$postIntervention == 1)
    points(tt,(yy), col = "#5D8AA8", pch = 19)
    points(tt[!intv],(yy[!intv]), col = "#FB9FA4", pch = 19)
    
    pop <-  stateInterventions$statePopulation[stateInterventions$stateAbbreviation ==   states[i]]
    
    bb <- posteriorMean$a[i]
    mm <- posteriorMean$b[i]
    ii <- min(which(intv)[1], 200, na.rm = T)
    
    curve(exp(bb + mm * x), 0, ii,col ="#FB9FA4", add = T)
    
    bb <- posteriorMean$a[i] + posteriorMean$g[i]
    mm <- posteriorMean$b[i] + posteriorMean$d[i]
    
    curve(exp(bb+mm*x), ii, 1000, col = "#5D8AA8", add = T)
    
  }
  if(!is.null(fileName)) {
    dev.off()
  }
}

#' @export
plotPropDeath <- function(allPropDead, plotT, endPlot, plotCol = plotCols[2]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qC <- t(apply(allPropDead, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  nomiss <- !is.na(rowSums(qC))
  plot(plotT, qC[, 2], col = NA, ylim = c(0, min(.12, max(qC[nomiss,3]))), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", xlab = "")
  
  polygon(c(plotT[nomiss], rev(plotT[nomiss])),
          c(qC[nomiss,1], rev(qC[nomiss,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT[nomiss], qC[nomiss,2], col = plotCol, lwd = 2)
  
  
  axisTicks <- seq(0,1,.02)
  
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = axisTicks, cex = .9, tick = F, hadj = .75)
  
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33", col   = "grey33", cex.axis = 1.25)
  
  mtext("Proportion of Resolved Cases Ending in Death", side=2, line=3.1, col="grey33", cex=1)
}

#' @export
plotDoubleTimeDeaths <- function(all2xDeaths, plotT, endPlot, plotCol = plotCols[2]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qC <- t(apply(all2xDeaths, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  nomiss <- !is.na(rowSums(qC))
  plot(plotT, qC[, 2], col = NA, ylim = c(0, max(qC[nomiss,3])), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", xlab = "")
  
  polygon(c(plotT[nomiss], rev(plotT[nomiss])),
          c(qC[nomiss,1], rev(qC[nomiss,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT[nomiss], qC[nomiss,2], col = plotCol, lwd = 2)
  
  
  axisTicks <- seq(0,max(qC, na.rm = T),100)
  
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33", col   = "grey33", cex.axis = 1.25)
  
  mtext("Cumulative Death Double Time (Days)", side=2, line=3.1, col="grey33", cex=1)
}

#' @export
plotDoubleTimeCases <- function(all2xCases, plotT, endPlot, plotCol = plotCols[2]) {
  par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
  qC <- t(apply(all2xCases, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  nomiss <- !is.na(rowSums(qC))
  plot(plotT, qC[, 2], col = NA, ylim = c(0, max(qC[nomiss,3])), xaxs = "i", xaxt = "n",
       yaxt = "n", ylab = "", xlab = "")
  
  polygon(c(plotT[nomiss], rev(plotT[nomiss])),
          c(qC[nomiss,1], rev(qC[nomiss,3])), col = paste0(plotCol,"30"), border = NA)
  lines(plotT[nomiss], qC[nomiss,2], col = plotCol, lwd = 2)
  
  
  axisTicks <- seq(0,max(qC, na.rm = T),100)
  
  axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
       labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  
  axisDays <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01",  "2020-12-01", "2021-03-01"))
  axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33", col   = "grey33", cex.axis = 1.25)
  
  mtext("Cumulative Case Double Time (Days)", side=2, line=3.1, col="grey33", cex=1)
}
