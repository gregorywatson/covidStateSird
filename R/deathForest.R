#' @export
deathForest <- function(stateCovidData, stateList, covariates, lagDays, fileOut = NULL) {
  
  velocLogCases <- velocLogDeaths <- data.frame()
  loc <- 0
  for(i in 1:length(stateList)) {
    loc <- loc + 1
     
    velocLoc <- velocitiesState(stateCovidData, stateList[i], stateInterventions, minCases = 0, endDate = endDate)
    covariatesLoc <- covariates[covariates$location == stateList[i],]
    population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == stateList[i]]
    velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases,  loc, covariatesLoc[,-1], population, row.names = NULL))
  }

  velocLogCasesList <- as.list(velocLogCases)
  velocLogCasesList$N <- nrow(velocLogCases)
  velocLogCasesList$nLoc <- length(unique(velocLogCasesList$loc))
  velocLogCasesList$p <- ncol(velocLogCases) - 3

  maxLag1 <- lagDays - 1
  deathModelData <- NULL
  for(i in unique(velocLogCases$loc)) {
    dDeaths <- diff(velocLogCases$deaths[velocLogCases$loc == i])[-(1:maxLag1)]
    dCases  <- diff(velocLogCases$u[velocLogCases$loc == i])
    
    laggedNewCases <- NULL
    for(k in 1:lagDays) {
      laggedNewCases <- cbind(laggedNewCases, lag(dCases,k-1)[-(1:maxLag1)])
    }
    colnames(laggedNewCases) <- paste0("dCase",1:lagDays)

    deathModelData <- rbind(deathModelData, cbind(velocLogCases[velocLogCases$loc ==   i,][-(1:(maxLag1+1)),], dDeaths, laggedNewCases))
  }
  
  randomForestDeathModel <- randomForest(dDeaths ~ ., data = deathModelData[,-(1:5)], na.action = na.roughfix, keep.forest=TRUE, keep.inbag=TRUE)
  
  if(!is.null(fileOut)) {
    save(randomForestDeathModel, deathModelData, file = fileOut)
  }
  return(randomForestDeathModel)
  
}
