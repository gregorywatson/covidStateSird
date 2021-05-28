#' @export
#' @importFrom randomForest randomForest na.roughfix
deathForest <- function(stateCovidData, stateList, covariates, lagDays, fileOut = NULL,
                        ntree = 500) {
  
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
    dDeaths <- diff(velocLogCases$deaths[velocLogCases$loc == i]) #[-(1:maxLag1)]
    dCases  <- diff(velocLogCases$u[velocLogCases$loc == i])
    
    laggedNewCases <- NULL
    laggedNewDeaths <- NULL
    for(k in 1:lagDays) {
      laggedNewCases <- cbind(laggedNewCases, lag(dCases,k-1)[-(1:maxLag1)])
      laggedNewDeaths <- cbind(laggedNewDeaths, lag(dDeaths,k-1)[-(1:maxLag1)])
    }
    colnames(laggedNewCases) <- paste0("dCase",1:lagDays)
    colnames(laggedNewDeaths) <- paste0("dDeaths",1:lagDays)
    
    deathModelData <- rbind(deathModelData, cbind(velocLogCases[velocLogCases$loc ==   i,][-(1:(maxLag1+1)),], dDeaths = dDeaths[-(1:maxLag1)], laggedNewCases, laggedNewDeaths))
  }
  
  randomForestDeathModel <- randomForest::randomForest(dDeaths ~ .,
                                                       data = deathModelData[,-c(1,3:5)], na.action = na.roughfix,
                                                       keep.forest=TRUE, keep.inbag=TRUE, ntree = ntree)
  
  if(!is.null(fileOut)) {
    save(randomForestDeathModel, deathModelData, file = fileOut)
  }
  return(randomForestDeathModel)
  
}


# from https://github.com/haozhestat/RFIntervals
# modified so doesn't retrain random forest
#' @export
RFOOBInterval <- function(rf,
                          x0,
                          alpha = 0.05,
                          symmetry = TRUE){
  
  ntest = nrow(x0)
  
  test_pred <- predict(rf, x0)
  oob_abs_error = sort(abs(rf$y - rf$predicted))
  oob_error = sort(rf$y - rf$predicted)
  
  upper_pred = rep(NA, ntest)
  lower_pred = rep(NA, ntest)
  
  ## symmetry = TRUE leads to the symmetric OOB Intervals
  ## symmetry = FALSE leads to the standard OOB Intervals
  if(symmetry){
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_abs_error,1-alpha)
      lower_pred[i] = test_pred[i] - quantile(oob_abs_error,1-alpha)
    }
  }
  else{
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_error, 1-alpha/2)
      lower_pred[i] = test_pred[i] + quantile(oob_error, alpha/2)
    }
    
  }
  
  return(list(pred = matrix(test_pred,ntest,1),
              lo = matrix(lower_pred,ntest,1),
              up = matrix(upper_pred,ntest,1),
              fit = matrix(rf$predicted,length(rf$predicted),1)))
  #   fit = matrix(predict(rf,x),n,1)))
}
