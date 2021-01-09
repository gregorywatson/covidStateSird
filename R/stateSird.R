#' @export
#' @importFrom dplyr %>% lag group_by mutate
#' @importFrom TTR runMean
stateSird <- function(stateAbbrev,
                      states,
                      covariates,
                      stateInterventions,
                      velocLogCases,
                      stateCovidData,
                      randomForestDeathModel,
                      posteriorSamples,
                      rfError = FALSE,
                      plots = TRUE,
                      makeRateTable = FALSE, 
                      n.t = 30,
                      lagDays = 21,
                      minCases = 100,
                      endDay = endDate,
                      endPlotDay = "2021-03-01") {

  statePop = stateInterventions$statePopulation[which(stateInterventions$stateAbbreviation == stateAbbrev)]

  stIndx <- which(states == stateAbbrev)

  nState <- sum(velocLogCases$loc == stIndx)
  u <- velocLogCases$u[which(velocLogCases$loc == stIndx)]
  y <- velocLogCases$y[which(velocLogCases$loc == stIndx)]
  
  n.day <- length(u)
  
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

  #n.t <- 201
  N.POST <- nrow(posteriorSamples$u)
  allStateFit <- NULL

  SAMPS <- 1:N.POST
  for(M in SAMPS) {

    nDay <-  length(days)
  
    gamm <- 1 / rnorm(1, 10, 1)
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
      steps[i, 4] <- steps[i-1, 4] + steps[i-1, 3] * gamm
      steps[i, 3] <- state[1, i] - steps[i, 4]
    }
    colnames(steps) <- c("times", "S", "I", "R")
  
   # init <- as.matrix(steps[nDay, 2:4, drop = F], nrow = 1)
   #  colnames(init) <- c("S", "I", "R")
    init <- as.numeric(steps[nDay, 2:4, drop = F])
    names(init) <- c("S", "I", "R")
    sirdParams <- list(
      phi = posteriorSamples$phi[M, stIndx],
      mu  = posteriorSamples$mu[M, stIndx],
      gamm = gamm,
      S0 = init[["S"]],
      Sminus1 = statePop - u[nState - 1],
      dSminus1 = -y[nState - 1] * u[nState - 1] * (statePop - u[nState - 1]) / init[["S"]],
      pop = statePop,
      tau = posteriorSamples$tau[M])
     
     stateFit <- data.frame(deSolve::dede(y = init, times = 0:(n.t-1), SIR, parms = sirdParams))
     # skip iteration if solution explodes
     if(sum(is.nan(stateFit$S))) next
     names(stateFit)[names(stateFit) == "time"] <- "times"
  
    stateFit$times <- days[nDay] + stateFit$times
  
  #  stateFit <- rbind(cbind(replicate = 1, steps[1:(nDay-1),]),  stateFit)
    stateFit <- cbind(replicate = 1, rbind(steps[1:(nDay-1),],  stateFit))
    
    # Estimate Deaths ---
    dCases  <- diff(statePop - stateFit$S)
  
    laggedNewCases <- NULL
     for(k in 1:lagDays) {
       laggedNewCases <- cbind(laggedNewCases, dplyr::lag(dCases,k-1)[-(1:(lagDays - 1))])
     }
     colnames(laggedNewCases) <- paste0("dCase",1:lagDays)
  
    x0 <- cbind(covariates[covariates$location == stateAbbrev,-1],
                loc = which(states == stateAbbrev),
                population = stateInterventions$statePopulation[
                  stateInterventions$stateAbbreviation == stateAbbrev],
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
    # ---
    
    
    stateFit$cases <- rowSums(stateFit[,c("I", "R", "D")])
    stateFit$doubleTimeCases <- stateFit$doubleTimeDeaths <- NA
                                         
    k <- 1
    l <- 1
    for(i in 1:nrow(stateFit)) {
      
      while((stateFit$cases[k] < 2 * stateFit$cases[i]) & (k <= nrow(stateFit))) {
        k <- k + 1
      }
      if((stateFit$cases[k] >= 2 * stateFit$cases[i]) & (k <= nrow(stateFit)))
        stateFit$doubleTimeCases[i] <- k
      
      while((stateFit$D[l] < 2 * stateFit$D[i]) & (l <= nrow(stateFit))) {
        l <- l + 1
      }
      if((stateFit$D[l] >= 2 * stateFit$D[i]) & (l <= nrow(stateFit)))
        stateFit$doubleTimeDeaths[i] <- l

    }
    # ---
    allStateFit <- rbind(allStateFit, stateFit)
  }

  nSampCompleted <- length(unique(allStateFit$replicate))
  nT <- nrow(steps) + n.t - 1

  allS <- matrix(allStateFit$S, nrow = nT,
                 ncol = nSampCompleted, byrow = F)
  allI <- matrix(allStateFit$I, nrow =  nT,
                 ncol = nSampCompleted, byrow = F)
  allR <- matrix(allStateFit$R, nrow =  nT,
                 ncol = nSampCompleted, byrow = F)
  allD <- matrix(allStateFit$D, nrow =  nT,
                 ncol = nSampCompleted, byrow = F)
  
  all2xCases <- matrix(allStateFit$doubleTimeCases, nrow =  nT,
                       ncol = nSampCompleted, byrow = F)
  all2xDeaths <- matrix(allStateFit$doubleTimeDeaths, nrow =  nT,
                        ncol = nSampCompleted, byrow = F)
  allStateFit$propDead <- allStateFit$D / (allStateFit$D + allStateFit$R)
  allPropDead <- matrix(allStateFit$propDead, nrow =  nT,
                        ncol = nSampCompleted, byrow = F)
    
  medianS <- as.matrix(apply(allS, 1, quantile, probs = c(.5),  na.rm = TRUE),
                       ncol = 1)
  medianI <- as.matrix(apply(allI, 1, quantile, probs = c(.5),  na.rm = TRUE),
                       ncol = 1)
  medianR <- as.matrix(apply(allR, 1, quantile, probs = c(.5),  na.rm = TRUE),
                       ncol = 1)
  medianD <- as.matrix(apply(allD, 1, quantile, probs = c(.5),  na.rm = TRUE),
                       ncol = 1)
  
 # nT <- length(stateFit$times)
  
  allStateFit <- as.data.frame(allStateFit %>% dplyr::group_by(replicate) %>%
                 dplyr::mutate(Rt = 10 * -c(diff(S), NA)   / I))
  
  allStateFit <- as.data.frame(allStateFit %>% dplyr::group_by(replicate) %>%
                 dplyr::mutate(Rt10 = c(TTR::runMean(10 *
                 - diff(S) / I[1:(nT-1)])[5:(nT-1)], rep(NA,5))))
  
  allRt <- matrix(allStateFit$Rt10, nrow = nT,
                  ncol = nSampCompleted, byrow = F)
  
  stateFit <- data.frame(cbind(c(steps[1:(nDay-1),1], days[nDay] + 0:(n.t-1)), medianS, medianI, medianR,
                               medianD))
  names(stateFit) <- c("times", "S", "I", "R", "D")
  stateFit[,1] <- as.Date(stateFit$times, origin = "1970-01-01")
  
  tableDays <- c(as.Date(endDate), as.Date(c("2021-02-01", "2021-03-01")))
  
  tRow <- which(stateFit$times %in% tableDays)

  if(makeRateTable) {
    rateTable <- data.frame(Day = c(as.Date(endDate),
                                    as.Date(c("2021-02-01", "2021-03-01"))),
               ActiveInfectionRate = 1e5 * (stateFit$I[tRow])/ statePop,
               DailyDeathRate = 1e5 * diff(stateFit$D)[tRow-1]/ statePop)
    write.csv(rateTable, file = paste0(outputPath, "/Data/", stateAbbrev,
                                       "_rateTable.csv"))
  }
  save(allStateFit, stateFit, sirdParams, file = paste0(outputPath, "/Data/",
       stateAbbrev, ".Rdata"))

  if(plots) {

    endPlot <- as.Date(endPlotDay)
    plotT <- stateFit$times[which(stateFit$times <= endPlot)]

    pdf(file = paste0(outputPath, "/Plots/",stateAbbrev,".pdf"),
        width = 24, height = 12)
    par(mfrow = c(2,4))
   
    plotCumulativeCases(allS, state, days, statePop,  plotT, endPlot)
    plotActiveCases(allI, state, plotT, endPlot)
    plotDailyDeaths(allD, allStateFit, stateFit, state, days, plotT, endPlot,
                    SAMPS, rfError, plotCol = plotCols[6])
    plotRt(allRt, state, days, plotT, endPlot)
    mtext(stateAbbrev, outer=TRUE,  cex=2, line=-4)
    
    plotDoubleTimeCases(all2xCases, plotT, endPlot, plotCols[7])
    plotDoubleTimeDeaths(all2xDeaths, plotT, endPlot, plotCols[8])
    plotPropDeath(allPropDead, plotT, endPlot, plotCols[3])
    
    dev.off()
  }
}


#' Multi-epoch SIRD
#'
#' Conducts simulations for a multi-epoch SIRD model
#'
#' @param n.rep The number of simulation replicates. Defaults to 10.
#' @param paramsList A list of parameter value lists. Each list entry must match \code{params} from \code{mcSeirScenario}.
#' @param init Initial values for \code{S}, \code{E}, \code{I}, \code{R} and if \code{func = seira}, \code{A}, the ambient source of infection.
#' @param timesList A list of time vectors.  Each list entry must match \code{times} from \code{mcSeirScenario} and be consecutive.
#' @param interventionList A list of intervention vectors. Each list entry must match \code{intervention} from \code{mcSeirScenario}.
#' @param nThreads The number of threads to use in multi-core computing.
#' @param func The function defining the derivatives in the ODE system (which define the model). Defaults to \code{seira}.
#' @importFrom tidyselect all_of
#' @export
multiEpochSird <- function(
  n.rep = 10,
  paramsList,
  timesList,
  interventionList,
  init,
  nThreads = 1,
  func = sird,
  modelOutput = if(is.list(init) & !is.data.frame(init)) names(init[[1]]) else names(init)
) {
  if (length(paramsList) != length(timesList) ||
      length(timesList) != length(interventionList)) {
    stop("Must specify the same number of epoches in `paramsList`, `timesList` and `interventionList`")
  }

  transformedTimesList <- transformEpochTimes(timesList)

  compartmentNames <- names(init)

  # First epoch
  result <- sirdModel(n.rep, paramsList[[1]], init, transformedTimesList[[1]],
                           interventionList[[1]], nThreads, func, model.output = modelOutput)

  last <- result
  if (length(paramsList) > 1) {
    for (e in 2:length(paramsList)) {
      df <- last$out %>% filter(.data$times == max(.data$times)) %>% arrange(replicate) %>% select(tidyselect::all_of(compartmentNames))
      init <- split(df, seq(nrow(df)))
      epoch <- sirdModel(n.rep, paramsList[[e]], init, transformedTimesList[[e]],
                              interventionList[[e]], nThreads, func, model.output = modelOutput)
      epoch$out$times <- epoch$out$times + min(timesList[[e]])
      last <- epoch
      result$out <- rbind(result$out, epoch$out[-1,])
    }
  }

  return(result)
}

transformEpochTimes <- function(timesList) {

  transformedTimeList <- list()
  startTime <- 0
  for (epoch in 1:length(timesList)) {
    times <- timesList[[epoch]]
    if (times[1] != startTime) {
      stop("Must specify consecutive times in `timesList`")
    }
    startTime <- times[length(times)]
    transformedTimeList[[epoch]] <- times - times[1]
  }

  return(transformedTimeList)
}

sirdModel <- function(
  n.rep = 10,
  params,
  init,
  times = seq(0, 365, by = 1),
  intervention = c(1, 1, 1),
  nThreads = 1,
  func = sird,
  model.output = if(is.list(init) & !is.data.frame(init)) names(init[[1]]) else names(init)
)
{
  n.t <- length(times)
  n.out <- n.t * n.rep
  
  # Set up output data frame -----------------------------------------------------
  if (is.list(init)) {
    template_init <- init[[1]]
  } else {
    template_init <- init
  }

  sirdOut <- data.frame(matrix(NA, ncol = 2 + length(model.output), nrow = n.out))

  for(j in 1:n.rep) {
    # Writing onto row [(j-1)*n.t+1] to row [j*n.t]:
    sirdOut[((j-1)*n.t + (1:n.t)),] <- execute.sird(j, n.t, params,
                                                        init, times, func)[, c("replicate", "time",   model.output)]
  }

  names(sirdOut) <- c("replicate", "times", model.output)

  sirdOut$intervention <- attr(intervention, "label")

  return(list(out = sirdOut))
}

execute.sird <- function(j, n.t, parameters, in_init, times, func) {

  if (is.list(in_init)) {
    init <- in_init[[j]]
    names <- names(init)
    init <- as.numeric(init)
    names(init) <- names
  } else {
    init <- in_init
  }

  result <- cbind(
    replicate = rep(j, n.t),
    as.data.frame(deSolve::ode(y = init, times = times, func = func, parms = parameters))
  )

  return(result)
}


#' SEIR Model
#'
#' This function defines the system of ODEs for a SEIR model. For use with deSolve::ode.
#'
#' @param times The times at which to evaluate the solutions.
#' @param state The current estimate of the variables in the ODE system.
#' @param parameters A vector of parameters. Must include \code{beta}, \code{rho}, \code{L} and \code{gamma}.
#' @export
seir <- function(times, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- - beta * S * I
    dE <-   beta * S * I - (E/L)
    dI <-                  (E/L) - gamma * I
    dR <-                          gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}

#' SIRD Model
#'
#' This function defines the system of ODEs for the SIRD model. For use with deSolve::ode.
#'
#' @param times The times at which to evaluate solutions.
#' @param state The current estimate of the variables in the ODE system.
#' @param parameters A vector of parameters. Must include parameters \code{t_s}, \code{a_s}, \code{b_s}, and \code{c_s}.
#' @export
sird <- function(times, state, parameters)
{
  with(as.list(c(state, parameters)), {
    du <- exp( (1/b_s) * exp(b_s * (times + t_s) + a_s) + c_s) * exp(b_s * (times + t_s) + a_s)
    dS <- - du
    dI <-   du - gamma * I
    dR <-   gamma * I
        return(list(c(dS, dI, dR)))
  })
}


#' @export
#' @importFrom deSolve lagvalue lagderiv
SIR <- function(time, state, parameters){
with(as.list(c(state,parameters)),{
  if(time < 2) {
    lagS <- Sminus1
    lagdS <- dSminus1
  }
  else {
    lagS  <- deSolve::lagvalue(time - 1, 1)
    lagdS <- deSolve::lagderiv(time - 1, 1)
  }
  
  u <- pop - S
  lag_u <- pop - lagS
  lagdu <- -lagdS * (S0/S)
   
  dS <- - ((lagdu / lag_u) ^ phi) * u * exp(mu + 1 / (2 * tau)) * (S/S0)
  dI <-   ((lagdu / lag_u) ^ phi) * u * exp(mu + 1 / (2 * tau)) * (S/S0)  - gamm * I
  dR <-                                                                     gamm * I
  return(list(c(dS,dI,dR)))})
}
