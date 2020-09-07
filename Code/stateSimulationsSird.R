library(TTR)
library(dplyr)
library(covidStateSird)

source(paste0(covidDir, "Code/McSeirScenario.R"))

# --- Train Death model --------------------------------------------------------------------
library(randomForest)

plotCols <- c("#AEC441", "#008875", "#0094D6", "#F2E000", "#C362A6", "#F78F1E")

stateCovidData <- read.csv(stateDataCSV)

covariates <- read.csv(paste0(covidDir, "Data/COVID-19 Location Covariates - analyticStates.csv"))
covariates <- covariates[, - which(names(covariates) == "stroke")]

covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] <-
covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] / 100000

set.seed(525600)
randomForestDeathModel <- deathForest(stateCovidData, states, covariates, 21, fileOut = paste0(outputPath, "/randomForestDeathModel.Rdata"))


stateInterventions <- read.csv(paste0(covidDir, "Data/StateInterventionDates.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)

load(paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))

set.seed(525600)
rm(list = "allStateFit")
stateSird("CA", covariates, stateInterventions, stateCovidData, randomForestDeathModel,
posteriorSamples, rfError = T)
load("/Users/gregw/Dropbox/Projects/covidStateSird/Output/States/2020-09-06/Data/CA.Rdata")
all.equal(allStateFit0, allStateFit)


foreach(i = 1:length(states))  %dopar% {
  stateRun(states[i], rfError = T)
}
