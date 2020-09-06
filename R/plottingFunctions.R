#' @export
plotVelocityFit <- function(posteriorMean, fileName = NULL) {
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
  
      bb <- pm$a[i]
      mm <- (pm$b[i])
      ii <- min(which(intv)[1], 200, na.rm = T)
  
      curve(exp(bb + mm * x), 0, ii,col ="#FB9FA4", add = T)
  
      bb <- pm$a[i] + pm$g[i]
      mm <- pm$b[i] + pm$d[i]
  
          curve(exp(bb+mm*x), ii, 1000, col = "#5D8AA8", add = T)
  
  }
  if(!is.null(fileName)) {
    dev.off()
  }
}

