#' @export
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

#' @export
riasJagsModel <- function(){
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

#' @export
constErrLogNorm <- function(x, t, u, beta, alpha) {
  err <- (u - exp((1/c(beta)) * exp(c(beta) * t + c(alpha)) + c(x))) ^2
  return(sum(err))
}

#' @export
caseModelConstant <- function(velocityPosterior, intervention = 1) {
  if(intervention == 0) {
    beta = velocityPosterior$BUGSoutput$sims.list$b
    alpha = velocityPosterior$BUGSoutput$sims.list$a
  } else if(intervention == 1) {
    beta = caseRiasFitMS$BUGSoutput$sims.list$b +
           caseRiasFitMS$BUGSoutput$sims.list$d
    alpha = caseRiasFitMS$BUGSoutput$sims.list$a +
            caseRiasFitMS$BUGSoutput$sims.list$g
  }
  
  constants <- foreach(k = 1:nrow(velocityPosterior$BUGSoutput$sims.matrix), .combine = 'rbind') %:%
     foreach(i = 1:velocLogCasesList$nLoc, .combine = 'c')  %dopar% {
       optim(6, constErrLogNorm, u = velocLogCases$u[(velocLogCases$loc == i) & (velocLogCases$postIntervention == intervention)],
            t = velocLogCases$t[(velocLogCases$loc == i) & (velocLogCases$postIntervention == intervention)],
         beta = beta[k,i],
         alpha = alpha[k,i],
       method = "SANN")$par
    }
  return(constants)
}
