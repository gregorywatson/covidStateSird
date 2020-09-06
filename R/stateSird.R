#' @export
stateSird <- function(stateAbbrev, covariates, stateInterventions, stateCovidData,
  randomForestDeathModel, posteriorSamples,
  rfError = FALSE, plots = TRUE, lagDays = 21, minCases = 100, endDay = endDate, endPlotDay = "2020-11-01") {

  statePop = stateInterventions$statePopulation[which(stateInterventions$stateAbbreviation == stateAbbrev)]

  stIndx <- which(states == stateAbbrev)

  stateLong <- stateCovidData[which(stateCovidData$state == stateAbbrev),]

  n.day <- nrow(stateLong)
  daysAll <- NULL
  for(i in 1:n.day) {
    daysAll[i] <-(paste(substr(stateLong$date[i],1,4),
                        substr(stateLong$date[i],5,6),
                        substr(stateLong$date[i],7,8), sep = "-") )
  }
  daysAll <- as.Date(daysAll)

  stateLong <- stateLong[order(daysAll),]
  daysAll   <- daysAll[order(daysAll)]

  stateAll <- rbind(stateLong$positive, stateLong$death, stateLong$hospitalizedCurrently)
  stateAll[which(is.na(stateAll))] <- 0
  
  startDay <- daysAll[which(stateAll[1,] >= minCases)][1]

  startDayIndx <- which(daysAll == startDay)

  if(!is.na(stateAll[2,startDayIndx])) {
    D0 <- stateAll[2,startDayIndx]
  } else D0 <- 0
  if(startDayIndx > lagDays) {
    R0 <- stateAll[1,startDayIndx - lagDays]
  } else R0 <- 0

  state <- stateAll[,which((daysAll >= startDay) & (daysAll <= endDay))]
  days  <- daysAll[which((daysAll >= startDay) & (daysAll <= endDay))]

  dataTimes <- days

  lastObs <- ncol(state)

  n.t <- 201
  N.POST <- nrow(posteriorSamples$a)
  allStateFit <- NULL

#  SAMPS <- sample(N.POST, 400) #1:N.POST
  SAMPS <- 1:N.POST
  for(M in SAMPS) {

    if(stateInterventions$interventionDate[which(stateInterventions$stateAbbreviation == stateAbbrev)]   == "") {
      a <- posteriorSamples$a[M, stIndx]
      b <- posteriorSamples$b[M, stIndx]
      c <- posteriorSamples$c_pre[M, stIndx]
    } else {
      a <- posteriorSamples$a[M, stIndx] + posteriorSamples$g[M, stIndx]
      b <- posteriorSamples$b[M, stIndx] + posteriorSamples$d[M, stIndx]
      c <- posteriorSamples$c_post[M, stIndx]
    }
    nDay <-  length(days)
  
    sirdParams <- c(
      scale = 1,
      t_s   = nDay,
      b_s   = b,
      a_s   = a,
      c_s   = c, 
      gamma = 1 / rnorm(1, 10, 1)
    )
  
    init0 <- matrix(c(S = statePop - stateAll[1,startDayIndx],
                       I = stateAll[1,startDayIndx] - (R0 + D0) ,
                       R = R0 + D0), nrow = 1)
    colnames(init0) <- c("S", "I", "R")
  
    steps <- data.frame(matrix(NA, nrow = nDay, ncol = 4))
    steps[,1] <- days
    steps[1, 2:4] <- init0
    steps[2:nDay, 2] <- init0[1] - state[1,2:nDay]
  
    for(i in 2:nDay) {
      steps[i, 2] <- statePop - state[1, i]
      steps[i, 4] <- steps[i-1, 4] + steps[i-1, 3] * sirdParams[["gamma"]]
      steps[i, 3] <- state[1, i] - steps[i, 4]
    }
    colnames(steps) <- c("times", "S", "I", "R")
  
    init <- as.matrix(steps[nDay, 2:4, drop = F], nrow = 1)
    colnames(init) <- c("S", "I", "R")
  
    stateFit <- try(mcMultiEpochSeirScenario(n.rep = 1,
                            params = list(sirdParams),
                            init = init[1,],
                            times = list(c(0:(n.t-1))),
                            intervention = list(c(1,1,1)),
                            func = sird,
                            modelOutput = c("S", "I", "R"))$out)
    if(inherits(stateFit, "try-error"))
    {
      # skip this iteration if there's an error (occasionally the ode solver doesn't converge)
      next
    }
  
    stateFit$times <- days[nDay] + stateFit$times
  
    stateFit <- rbind(cbind(replicate = 1, steps[1:(nDay-1),]),  stateFit)
  
    # Estimate Deaths ---
    dCases  <- diff(statePop - stateFit$S)
  
    laggedNewCases <- NULL
     for(k in 1:lagDays) {
       laggedNewCases <- cbind(laggedNewCases, lag(dCases,k-1)[-(1:(lagDays - 1))])
     }
     colnames(laggedNewCases) <- paste0("dCase",1:lagDays)
  
    x0 <- cbind(covariates[covariates$location == stateAbbrev,-1],
                loc = which(states == stateAbbrev),
                population = stateInterventions$statePopulation[stateInterventions$stateAbbreviation ==   stateAbbrev],
                laggedNewCases, row.names = NULL)
  
    if(rfError == TRUE) {
      rfInt <- RFOOBInterval(randomForestDeathModel, x0 = x0)
      deltaDeathPred <- rfInt$pred
      rfLo <- rfInt$lo
      rfUp <- rfInt$up
    } else {
      deltaDeathPred <- predict(randomForestDeathModel, newdata = x0)
      rfLo <-deltaDeathPred
      rfUp <-deltaDeathPred
    }
    newCaseSum <- rowSums(x0[,paste0("dCase", 1:lagDays)])
    pctCap <- c(rep(.15, 30), rep(.08, length(newCaseSum) - 30))
  
    deltaDeathPredCap <- pmin(.07 * pctCap * newCaseSum, deltaDeathPred)
    rfLoCap <- pmax(0, rfLo)
    rfLoCap <- pmin(.01 * pctCap * newCaseSum, rfLoCap)
    rfUpCap <- pmax(0, rfUp)
    rfUpCap <- pmin(.15 * pctCap * newCaseSum, rfUpCap)
  
    deltaDeath <-  c(rep(0, lagDays), deltaDeathPredCap)
  
    if(ncol(state) >= lagDays) {
      deltaDeath[2:lagDays] <- diff(state[2,1:lagDays])
    } else {
      deltaDeath[2:lagDays] <- 0
    }
    deltaDeath[is.na(deltaDeath)] <- 0
  
    stateFit$D <- rep(0, length(stateFit$S))
    stateFit$D <- cumsum(deltaDeath)
  
    stateFit$R <- stateFit$R - stateFit$D
  
    stateFit$replicate <- M
  
    stateFit$deathPred <- deltaDeath
    if(rfError) {
      stateFit$deathLo <- c(deltaDeath[1:lagDays], rfLoCap)
      stateFit$deathUp <- c(deltaDeath[1:lagDays], rfUpCap)
    }
    allStateFit <- rbind(allStateFit, stateFit)
  }


  allS <- matrix(allStateFit$S, nrow = length(stateFit$times), ncol =length(SAMPS), byrow = F)
  allI <- matrix(allStateFit$I, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
  allR <- matrix(allStateFit$R, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
  allD <- matrix(allStateFit$D, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
  
  medianS <- as.matrix(apply(allS, 1, quantile, probs = c(.5),  na.rm = TRUE), ncol = 1)
  medianI <- as.matrix(apply(allI, 1, quantile, probs = c(.5),  na.rm = TRUE), ncol = 1)
  medianR <- as.matrix(apply(allR, 1, quantile, probs = c(.5),  na.rm = TRUE), ncol = 1)
  medianD <- as.matrix(apply(allD, 1, quantile, probs = c(.5),  na.rm = TRUE), ncol = 1)
  
  
  nT <- length(stateFit$times)
  
  allStateFit <- as.data.frame(allStateFit %>% group_by(replicate) %>% mutate(Rt = 10 * -c(diff(S), NA)   / I))
  
  allStateFit <- as.data.frame(allStateFit %>% group_by(replicate) %>% mutate(Rt10 = c(runMean(10 *   -diff(S) / I[1:(nT-1)])[5:(nT-1)], rep(NA,5))))
  
  allRt <- matrix(allStateFit$Rt10, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
  
  stateFit <- data.frame(cbind(stateFit$times, medianS, medianI, medianR, medianD))
  names(stateFit) <- c("times", "S", "I", "R", "D")
  stateFit[,1] <- as.Date(stateFit$times, origin = "1970-01-01")
  
  tableDays <- c(as.Date(endDate), as.Date(c("2020-10-01", "2020-11-01")))
  
  tRow <- which(stateFit$times %in% tableDays)
  
  rateTable <- data.frame(Day = c(as.Date(endDate), as.Date(c("2020-10-01", "2020-11-01"))),
             ActiveInfectionRate = 1e5 * (stateFit$I[tRow])/ statePop,
             DailyDeathRate = 1e5 * diff(stateFit$D)[tRow-1]/ statePop)
  write.csv(rateTable, file = paste0(outputPath, "/Data/", stateAbbrev, "_rateTable.csv"))
  
  save(allStateFit, stateFit, pTuned, file = paste0(outputPath, "/Data/", stateAbbrev, ".Rdata"))

  if(plots) {

    endPlot <- as.Date(endPlotDay)
    plotT <- stateFit$times[which(stateFit$times <= endPlot)]
    
    pdf(file = paste0(outputPath, "/Plots/",stateAbbrev,".pdf"), width = 24,   height = 6)
    par(mfrow = c(1,4))
    
    plotCumulativeCases(allS, state, statePop,  plotT, endPlot)
    plotActiveCases(allI, state, plotT, endPlot)
    plotDailyDeaths(allD, allStateFit, stateFit, state, plotT, endPlot, SAMPS, rfError, plotCol = plotCols[6])
    plotRt(allRt, state, plotT, endPlot)

    dev.off()
  }
}
