library(R2jags)
library(dplyr)
library(randomForest)
library(doParallel)

devtools::load_all()

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

plotVelocityFit(caseRiasFitMS$BUGSoutput$mean,
                fileName = paste0(outputPath, "/Plots/velocityModelFit.pdf"))

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

