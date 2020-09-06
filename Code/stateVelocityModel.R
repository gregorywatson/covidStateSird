library(R2jags)
library(dplyr)
library(randomForest)
library(doParallel)

# the last day of data to use
endDate <- "2020-08-30"

covidDir <- "/Users/gregw/Dropbox/Projects/covidStateSird/"

registerDoParallel(cores=5)

outputPath <- file.path(covidDir, "Output/States", Sys.Date())
dir.create(outputPath)
dir.create(file.path(outputPath, "/Tables/"))
dir.create(file.path(outputPath, "/Plots/"))
dir.create(file.path(outputPath, "/Data/"))

# ===============================================================

velocitiesState <- function(statesLong, stateAbbrev, stateInterventions = NULL, minCases = 100, endDate = "2099-01-01") {

  if(!is.null(stateInterventions)) {
    intervDate <- stateInterventions$interventionDate[stateInterventions$stateAbbreviation == stateAbbrev]
    if(intervDate == "") intervDate <- "2099-01-01"
  }
  stateLong <- statesLong[which(statesLong$state == stateAbbrev),]
  n.day <- nrow(stateLong)
  daysAll <- NULL
  for(i in 1:n.day) {
    daysAll[i] <-(paste(substr(stateLong$date[i],1,4),
                        substr(stateLong$date[i],5,6),
                        substr(stateLong$date[i],7,8), sep = "-") )
  }
  daysAll <- as.Date(daysAll)

  stateLong <- stateLong[order(daysAll, decreasing = F),]
  daysAll   <- daysAll[order(daysAll, decreasing = F)]
        
  stateAll <- rbind(stateLong$positive, stateLong$death, stateLong$hospitalizedCurrently)
  state <- stateAll[,which(stateAll[1,] >= minCases)]
  days  <- as.Date(daysAll[which(stateAll[1,] >= minCases)])
  
  postIntervention <- 1 * (days > intervDate)
  
  y <- log(state[1, !is.na(state[1,]) & (days <= endDate)])
  x <- which(!is.na(state[1,]) & (days <= endDate))
  splineCases <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  derivCases <- predict(splineCases, deriv = 1)
  
  y <- log(state[2, !is.na(state[2,]) & (days <= endDate)])
  x <- which(!is.na(state[2,]) & (days <= endDate))
  splineDeaths <- smooth.spline(x = x[y > - Inf], y = y[y > -Inf])
  derivDeaths <- predict(splineDeaths, deriv = 1)

  return(list(cases = data.frame(y = derivCases$y,
                                 t = derivCases$x,
                                 u = state[1,derivCases$x],
                                 postIntervention = postIntervention[derivCases$x],
                                 deaths = state[2, derivCases$x]),
              deaths = data.frame(y = derivDeaths$y,
                                  t = derivDeaths$x,
                                  u = state[2,derivDeaths$x],
                                  hosps = state[3,derivDeaths$x],
                                 postIntervention = postIntervention[derivDeaths$x]),
              days = days,
              intervDate = intervDate))
}


stateInterventions <- read.csv(paste0(covidDir, "Data/StateInterventionDates.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)

stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"

statesLong <- read.csv(stateDataCSV)

covariates <- read.csv(paste0(covidDir, "Data/COVID-19 Location Covariates - analyticStates.csv"))
# covariates <- covariates[complete.cases(covariates),]
covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] <-
  covariates[, c("heartDisease", "lungCancer", "diabetes", "stroke", "copd")] / 100000

velocLogCases <- velocLogDeaths <- data.frame()
loc <- 0

states <- c("NY", "CA", "LA", "NJ", "MI", "IL", "GA", "FL", "TX", "PA", "AK", "AZ", "CO", "IN", "MN", "WA", "MO", "NV", "OH", "VA", "WI", "OR", "AL", "AR", "DE", "HI", "ID", "IA", "KS", "KY", "ME", "MS", "MT", "NE", "NM", "NC", "ND", "OK", "RI", "SC", "SD", "TN", "UT", "WV", "CT", "MD", "MA",
    "NH", "VT", "WY")
    

for(i in 1:length(states)) {
  loc <- loc + 1

  if(states[i] == "SD")   minCase <- 102
  else minCase <- 100
  
  velocLoc <- velocitiesState(statesLong, states[i], stateInterventions, minCases = minCase, endDate = endDate)
  covariatesLoc <- covariates[covariates$location == states[i],]
  population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == states[i]]
  
  velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases,  loc, covariatesLoc[,-1], population))
}

velocLogCasesList <- as.list(velocLogCases)
velocLogCasesList$N <- nrow(velocLogCases)
velocLogCasesList$nLoc <- length(unique(velocLogCasesList$loc))
velocLogCasesList$p <- ncol(velocLogCases) - 3

velocLogCasesListMS <- velocLogCasesList
velocLogCasesListMS$y[velocLogCasesList$y <= 0] <- NA
velocLogCasesListMS$p <- 12

riasJagsMS <- function(){
  # Likelihood:
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau[i])

    mu[i] <-  a[loc[i]] + b[loc[i]] * t[i] + (g[loc[i]] + d[loc[i]] * t[i]) * postIntervention[i]
    tau[i] <- exp(alpha[loc[i]] + beta[loc[i]] * t[i])
  }
  
  # Priors:
  for(j in 1:nLoc) {
    a[j] ~ dnorm(mu_a, 1) # random intercept for location
    b[j] ~ dnorm(mu_b, 10) # random slope for location
    g[j] ~ dnorm(mu_g, 1) # random effect for post-intervention
    d[j] ~ dnorm(mu_d, 10) # random slope for post-intervention
    alpha[j] ~ dnorm(mu_alpha, 1)
       beta[j]  ~ dnorm(mu_beta, 10)
  }
  mu_a  ~ dnorm(-1, 1) # random intercept mean
  mu_b  ~ dnorm(0, 10) # random slope mean
  mu_g  ~ dnorm(0, 1) # random effect mean
  mu_d  ~ dnorm(-.05, 10) # random slope mean
  mu_alpha ~ dnorm(0,1)
  mu_beta  ~ dnorm(0, 10)
}

 params <- c("mu_a", "mu_b", "sigma", "a", "b", "g", "d", "mu_g", "mu_d", "mu", "alpha", "beta", "mu_alpha", "mu_beta")
 
 start <- Sys.time()

caseRiasFitMS <- jags.parallel(data = velocLogCasesListMS, inits = NULL, parameters.to.save = params,
  model.file = riasJagsMS, n.chains = 3, n.iter = 200, n.burnin = 10, n.thin = 20, DIC = F)

 timeElapsed <- (Sys.time() - start)

pm <- caseRiasFitMS$BUGSoutput$mean

pdf(file = paste0(outputPath, "/Plots/velocityModelFit.pdf"), width = 6, height = 6)
for(i in 1:length(states)) {
    testV <- velocitiesState(statesLong, states[i], stateInterventions, minCases = 100)
    tt <- testV$cases$t
    y <- testV$cases$y
    y[which(y <= 0)] <- NA
    yy <-(y)
    
  plot(tt,(yy), col = "white", pch = 19, type = "o", lwd = .5, cex = .7,
       main = states[i], ylim = range((yy), na.rm = T), ylab = "d/dt Log Cumulative Cases", xlab = "Days Since 100+ Cases")
    intv <- (testV$cases$postIntervention == 1)
    points(tt,(yy), col = "#5D8AA8", pch = 19)
    points(tt[!intv],(yy[!intv]), col = "#FB9FA4", pch = 19)
    
    pop <-  stateInterventions$statePopulation[stateInterventions$stateAbbreviation == states[i]]
    
    bb <- pm$a[i]
    mm <- (pm$b[i])
    ii <- min(which(intv)[1], 200, na.rm = T)

    curve(exp(bb + mm * x), 0, ii,col ="#FB9FA4", add = T)
    
    bb <- pm$a[i] + pm$g[i]
    mm <- pm$b[i] + pm$d[i]

        curve(exp(bb+mm*x), ii, 1000, col = "#5D8AA8", add = T)

}
dev.off()

constErrLogNorm <- function(x, t, u, beta, alpha) {
  err <- (u - exp((1/c(beta)) * exp(c(beta) * t + c(alpha)) + c(x))) ^2
  return(sum(err))
}

 c_S_pre <- foreach(k = 1:nrow(caseRiasFitMS$BUGSoutput$sims.matrix), .combine = 'rbind') %:%
   foreach(i = 1:velocLogCasesList$nLoc, .combine = 'c')  %dopar% {
     optim(6, constErrLogNorm, u = velocLogCases$u[(velocLogCases$loc == i) & (velocLogCases$postIntervention == 0)],
          t = velocLogCases$t[(velocLogCases$loc == i) & (velocLogCases$postIntervention == 0)],
       beta = caseRiasFitMS$BUGSoutput$sims.list$b[k,i],
       alpha = caseRiasFitMS$BUGSoutput$sims.list$a[k,i],
     method = "SANN")$par
}
 
c_S_post <- foreach(k = 1:nrow(caseRiasFitMS$BUGSoutput$sims.matrix), .combine = 'rbind') %:%
   foreach(i = 1:velocLogCasesList$nLoc, .combine = 'c')  %dopar% {
     optim(6, constErrLogNorm, u = velocLogCases$u[velocLogCases$loc == i],
       t = velocLogCases$t[velocLogCases$loc == i],
       beta = caseRiasFitMS$BUGSoutput$sims.list$b[k,i] + caseRiasFitMS$BUGSoutput$sims.list$d[k,i],
       alpha = caseRiasFitMS$BUGSoutput$sims.list$a[k,i] + caseRiasFitMS$BUGSoutput$sims.list$g[k,i],
       method = "SANN") $par
}
rownames(c_S_pre) <- NULL
rownames(c_S_post) <- NULL

posteriorSamples <- caseRiasFitMS$BUGSoutput$sims.list

posteriorSamples[["c_pre"]] <- c_S_pre
posteriorSamples[["c_post"]] <- c_S_post

save(posteriorSamples,
  file = paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))

# run state SIRD models
source(paste0(covidDir, "Code/stateSimulationsSird.R"))

