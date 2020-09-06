


library(TTR)
library(dplyr)

source(paste0(covidDir, "Code/McSeirScenario.R"))
source(paste0(covidDir, "Code/RFOOBInterval.R"))

# --- Train Death model --------------------------------------------------------------------
library(randomForest)
library(dplyr) # critical for lag function

plotCols <- c("#AEC441", "#008875", "#0094D6", "#F2E000", "#C362A6", "#F78F1E")

covariates <- read.csv(paste0(covidDir, "Data/COVID-19 Location Covariates - analyticStates.csv"))

covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] <-
  covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] / 100000
  
statesRF <- c("NY", "CA", "LA", "NJ", "MI", "IL", "GA", "FL", "TX", "PA", "AK", "AZ", "CO", "IN", "MN", "WA", "MO", "NV", "OH", "VA", "WI", "OR", "AL", "AR", "DE", "HI", "ID", "IA", "KS", "KY", "ME", "MS", "MT", "NE", "NM", "NC", "ND", "OK", "RI", "SC", "SD", "TN", "UT", "WV", "CT", "MD", "MA",
    "NH", "VT", "WY")


velocLogCases <- velocLogDeaths <- data.frame()
loc <- 0
for(i in 1:length(statesRF)) {
  loc <- loc + 1
   
  velocLoc <- velocitiesState(statesLong, statesRF[i], stateInterventions, minCases = 0, endDate = endDate)
  covariatesLoc <- covariates[covariates$location == statesRF[i],]
  population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == statesRF[i]]
  velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases,  loc, covariatesLoc[,-1], population))
}

velocLogCasesList <- as.list(velocLogCases)
velocLogCasesList$N <- nrow(velocLogCases)
velocLogCasesList$nLoc <- length(unique(velocLogCasesList$loc))
velocLogCasesList$p <- ncol(velocLogCases) - 3

lagDays <- 21
maxLag1 <- lagDays - 1
deathModel <- NULL
for(i in unique(velocLogCases$loc)) {
  dDeaths <- diff(velocLogCases$deaths[velocLogCases$loc == i])[-(1:maxLag1)]
  dCases  <- diff(velocLogCases$u[velocLogCases$loc == i])
  
  laggedNewCases <- NULL
  for(k in 1:lagDays) {
    laggedNewCases <- cbind(laggedNewCases, lag(dCases,k-1)[-(1:maxLag1)])
  }
  colnames(laggedNewCases) <- paste0("dCase",1:lagDays)

  deathModel <- rbind(deathModel, cbind(velocLogCases[velocLogCases$loc ==   i,][-(1:(maxLag1+1)),], dDeaths, laggedNewCases))
}
  
# remove strokes as covariate due to missing values
deathModel <- deathModel[, - which(names(deathModel) == "stroke")]
# deathModel <- deathModel[, - which(names(deathModel) == "loc")]
set.seed(525600)
randomForestDeathModel <- randomForest(dDeaths ~ ., data = deathModel[,-(1:5)], na.action = na.roughfix, keep.forest=TRUE, keep.inbag=TRUE)

save(randomForestDeathModel, deathModel, file = paste0(outputPath, "/randomForestDeathModel.Rdata"))

# ----------------------------------------------------------------------------------------------

states <- c("NY", "CA", "LA", "NJ", "MI", "IL", "GA", "FL", "TX", "PA", "AK", "AZ", "CO", "IN", "MN", "WA", "MO", "NV", "OH", "VA", "WI", "OR", "AL", "AR", "DE", "HI", "ID", "IA", "KS", "KY", "ME", "MS", "MT", "NE", "NM", "NC", "ND", "OK", "RI", "SC", "SD", "TN", "UT", "WV", "CT", "MD", "MA", "NH", "VT", "WY")


SAMPS <- sample(2850, 400)
stateRun <- function(stateAbbrev, rfError = FALSE, plots = TRUE) {

  minCases <- 100
  endDay <- endDate

  set.seed(525600)

  stateInterventions <- read.csv(paste0(covidDir, "Data/StateInterventionDates.csv"),
                                 stringsAsFactors = FALSE,
                                 header = TRUE)

  statePop = stateInterventions$statePopulation[which(stateInterventions$stateAbbreviation == stateAbbrev)]

  stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"

  load(paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))
  stIndx <- which(states == stateAbbrev)

  load(file = paste0(outputPath, "/randomForestDeathModel.Rdata"))
 
  covariates <- read.csv(paste0(covidDir, "/Data/COVID-19 Location Covariates - analyticStates.csv"))
  covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] <-
    covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] / 100000
  

  stateLong <- read.csv(stateDataCSV)
  stateLong <- stateLong[which(stateLong$state == stateAbbrev),]

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
  
    pSir2 <- c(
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
      steps[i, 4] <- steps[i-1, 4] + steps[i-1, 3] * pSir2[["gamma"]]
      steps[i, 3] <- state[1, i] - steps[i, 4]
    }
    colnames(steps) <- c("times", "S", "I", "R")
  
    init <- as.matrix(steps[nDay, 2:4, drop = F], nrow = 1)
    colnames(init) <- c("S", "I", "R")
  
  
    pTuned <- pSir2
  
    stateFit <- try(mcMultiEpochSeirScenario(n.rep = 1,
                            params = list(pTuned),
                            init = init[1,],
                            times = list(c(0:(n.t-1))),
                            intervention = list(c(1,1,1)),
                            func = sir7,
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
       laggedNewCases <- cbind(laggedNewCases, lag(dCases,k-1)[-(1:maxLag1)])
     }
     colnames(laggedNewCases) <- paste0("dCase",1:lagDays)
  
    x0 <- cbind(covariates[covariates$location == stateAbbrev,-1],
                loc = which(statesRF == stateAbbrev),
                population = stateInterventions$statePopulation[stateInterventions$stateAbbreviation ==   stateAbbrev],
                laggedNewCases, row.names = NULL)
  
    x0 <- x0[, - which(names(x0) == "stroke")]
  
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
  
  
  tableDays <- c(as.Date(endDate), as.Date(c("2020-09-01", "2020-10-01")))
  
  tRow <- which(stateFit$times %in% tableDays)
  
  rateTable <- data.frame(Day = c(as.Date(endDate), as.Date(c("2020-09-01", "2020-10-01"))),
             ActiveInfectionRate = 1e5 * (stateFit$I[tRow])/ statePop,
             DailyDeathRate = 1e5 * diff(stateFit$D)[tRow-1]/ statePop)
  write.csv(rateTable, file = paste0(outputPath, "/Data/", stateAbbrev, "_rateTable.csv"))
  
  save(allStateFit, stateFit, pTuned, file = paste0(outputPath, "/Data/", stateAbbrev, ".Rdata"))


  if(plots) {
  
    # --- cumulative cases posterior -------------------------------------------------
  
    plotCols <- c("#AEC441", "#008875", "#0094D6", "#F2E000", "#C362A6", "#F78F1E")
  
    endPlot <- as.Date("2020-11-01")
    plotT <- stateFit$times[which(stateFit$times <= endPlot)]
  
  pdf(file = paste0(outputPath, "/Plots/",stateAbbrev,".pdf"), width = 24,   height = 6)
  par(mfrow = c(1,4))

    par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
    qC <- t(apply(statePop - allS, 1, quantile, probs = c(.05, .5,  .95),  na.rm =     TRUE))[1:length(plotT),]
  
    plot(plotT, qC[,2], col = NA, ylim = c(0, max(qC[,3])), xaxs = "i", xaxt = "n",
         yaxt = "n", ylab = "", xlab = "")
  
    polygon(c(plotT, rev(plotT)),
    c(qC[,1], rev(qC[,3])), col = paste0(plotCols[2],"30"), border = NA)
    lines(plotT, qC[,2], col = plotCols[2], lwd = 2)
    points(days, state[1,], pch = 19, col = "grey33")
  
    if(max(qC) > 500000) {
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

    axisDays <- as.Date(c("2020-03-01", "2020-05-01", "2020-07-01",  "2020-09-01", "2020-11-01"))
   axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33", col   = "grey33", cex.axis = 1.25)

     mtext("Cumulative Confirmed Cases", side=2, line=3.1, col="grey33", cex=1)
  
    # --- active cases posterior ----------
  
    par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
    qI <- t(apply(allI, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
    plot(plotT, qI[,2], col = NA, ylim = c(0, max(qI[,3]) * 1.04), xaxs = "i", xaxt = "n",
          yaxt = "n", ylab = "", yaxs = "i", xlab = "")
  
    polygon(c(plotT, rev(plotT)),
    c(qI[,1], rev(qI[,3])), col = paste0(plotCols[5],"30"), border = NA)
    lines(plotT, qI[,2], col = plotCols[5], lwd = 2)
    #points(days, state[1,], pch = 19, col = "grey33")
  
    if(max(qI) > 100000) {
      axisTicks <- seq(0,500000,25000)
    } else if (max(qI) > 10000) {
      axisTicks <- seq(0,500000,5000)
    } else if (max(qI) > 1000) {
      axisTicks <- seq(0,500000,500)
    } else if (max(qI) > 100) {
      axisTicks <- seq(0,500000,50)
    } else {
      axisTicks <- seq(0,500000,5)
    }
      axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
     labels = formatC(axisTicks, format = "d", big.mark = ","), cex = .9, tick = F, hadj = .75)
  
   axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
     mtext("Actively Infected Confirmed Cases", side=2, line=3.1, col="grey33", cex=1)
  
    # --- daily deaths posterior --------------------
    par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
    allDD <- apply(allD, 2, diff)
    qDD <- t(apply(allDD, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
    if(rfError) {
      wideDeathLo <- matrix(allStateFit$deathLo, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
      qDD[,1] <- t(apply(wideDeathLo, 1, quantile, probs = c(.05,.5),  na.rm = TRUE))[1:length(plotT),1]
      wideDeathUp <- matrix(allStateFit$deathUp, nrow = length(stateFit$times), ncol = length(SAMPS), byrow = F)
      qDD[,3] <- t(apply(wideDeathUp, 1, quantile, probs = c(.95,.5),  na.rm = TRUE))[1:length(plotT),1]
    }
  
    plot(plotT, qDD[,2], col = NA, ylim = c(0, max(qDD[,3]) * 1.04), xaxs = "i", xaxt = "n",
          yaxt = "n", ylab = "", yaxs = "i", xlab = "")
  
    polygon(c(plotT, rev(plotT)),
    c(qDD[,1], rev(qDD[,3])), col = paste0(plotCols[6],"30"), border = NA)
    lines(plotT, qDD[,2], col = plotCols[6], lwd = 2)
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
  
    axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
     mtext("Daily Deaths", side=2, line=3.1, col="grey33", cex=1)
  
    # Rt posterior plot -----
    par(mai = c(.8,.8,1,.4), mgp = c(3,.75,0))
  
    qRt <- t(apply(allRt, 1, quantile, probs = c(.05, .5,  .95),  na.rm = TRUE))[1:length(plotT),]
  
    plot(plotT, qRt[,2], col = NA, ylim = c(0, max(qRt[,3], na.rm = T)), xaxs = "i",
         xaxt = "n", yaxt = "n", ylab = "", xlab = "")
    abline(h = 1, col = "grey66", lty = 2)
    polygon(c(plotT, rev(plotT)),
    c(qRt[,1], rev(qRt[,3])), col = paste0(plotCols[1],"30"), border = NA)
    lines(plotT, qRt[,2], col = plotCols[1], lwd = 2)
    points(days, state[1,], pch = 19, col = "grey33")
  
    axisTicks <- seq(0,4,.5)
    axis(2, at = axisTicks, col = "grey33", las = 2, col.axis = "grey33",
    labels = axisTicks, cex = .9, tick = F, hadj = .75)
  
    mtext("Estimated Rt 10-day Moving Average", side=2, line=3.1, col="grey33", cex=1)
  
    axis(1, at = axisDays, labels = gsub(" 0", " ", format(axisDays, "%B %d")), col.axis = "grey33",   col   = "grey33", cex.axis = 1.25)
  
    mtext(stateAbbrev, outer=TRUE,  cex=2, line=-4)
    dev.off()
  }
}


foreach(i = 1:length(states))  %dopar% {
  stateRun(states[i], rfError = T)
}

states2 <- c("MA", "MD", "TN", "SD", "ND", "NC", "ME", "HI", "WI", "WA")

for(i in 1:length(states2)) {
  stateRun(states2[i], rfError = T)
}
