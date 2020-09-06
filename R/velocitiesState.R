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
