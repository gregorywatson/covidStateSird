
# ------ RF importance score ------
par(mai = c(1.5,0,0,0))
nVar <- length(rownames(randomForestDeathModel$importance))
barplot(t(randomForestDeathModel$importance), col = plotCols[2], border = NA,
ylab = "", yaxt = "n", xaxt = "n")


varnames <- c("Location", "Male", "Age 0-24", "Age 25-34", "Age 35-44", "Age 45-54", "Age 55-64",
  "Age 65+", "Heart Disease", "Lung Cancer", "Diabetes", "COPD", "Urbanization", "Population", paste0("New Cases on t-",1:lagDays))


axis(1, at = 1.2 * (1:nVar) -.5 , labels = varnames, las = 3, tick = F,
   cex.axis = 1, line = -1)

mtext("Normalized Variable Importance", side=2, line=-2, col="grey33")

# --- log scale
pdf(file = paste0(outputPath, "/Plots/Death_Model_Log_Importance.pdf"), width = 12, height = 7)
par(mai = c(1.5,0,0,0))
nVar <- length(rownames(randomForestDeathModel$importance))

barplot(c(t(log(as.numeric(randomForestDeathModel$importance)))), border = NA,
ylab = "", yaxt = "n", xaxt = "n",
col = c(paste0(plotCols[5],"40"), paste0(plotCols[3],"40"), rep(paste0(plotCols[2],"40"),6), rep(paste0(plotCols[1],"40"),4), rep(paste0(plotCols[4],"40"),2), rep(paste0(plotCols[6],"40"),lagDays))
)

varnames <- c("Location", "Male", "Age 0-24", "Age 25-34", "Age 35-44", "Age 45-54", "Age 55-64",
 "Age 65+", "Heart Disease", "Lung Cancer", "Diabetes", "COPD", "Urbanization", "Population", paste0("New Cases on t-",1:lagDays))
axis(1, at = 1.2 * (1:nVar) -.5 , labels = varnames, las = 3, tick = F,
cex.axis = 1, line = -.75, col.axis = "grey33")

mtext("Normalized Variable Importance (Log Scale)", side=2, line=-1.5, col="grey33")
dev.off()





# --- prediction error ---

sirdPredError <- function(testDays, dataDir) {
  caseError <- deathError <- cumulCase <- predCumulCase <-
    cumulDeath <- predCumulDeath <- newCase <- newDeath <-
    predNewCase <- predNewDeath <- matrix(NA, nrow = 50, ncol = length(testDays))

  names(caseError) <- names(deathError) <- c(as.character(testDays))

  maseDenom <- rep(NA< 50)

  for(j in 1:length(states)) {
    stateAbbrev <- states[j]

    statePop <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation ==   stateAbbrev]

    load(paste0(dataDir, stateAbbrev, ".Rdata"))

    stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"
  
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
  
    predDays <- which(stateFit$times %in% testDays)
    obsDays <- which(daysAll %in% testDays)
  
    newCase[j,]  <- c(NA,diff(stateLong$positive))[obsDays]
    newDeath[j, ]<- c(NA,diff(stateLong$death))[obsDays]
  
    predNewCase[j,] <- c(NA,diff(statePop - stateFit$S))[predDays]
    predNewDeath[j,]  <- c(NA, diff(stateFit$D))[predDays]
  
    cumulCase[j,] <- stateLong$positive[obsDays]
    predCumulCase[j,] <- (statePop - stateFit$S)[predDays]
  
    cumulDeath[j,] <- stateLong$death[obsDays]
    predCumulDeath[j,] <- stateFit$D[predDays]
  
    maseDenomCase <- mean(abs(diff(newCase)), na.rm = T)
    maseDenomDeath <- mean(abs(diff(newDeath)), na.rm = T)
  }
 
  caseError    <- abs(newCase - predNewCase) / maseDenomCase
  deathError   <- abs(newDeath - predNewDeath) / maseDenomDeath

  nTestDays <- length(testDays)

  meanCaseError <- rep(NA, nTestDays)
  meanDeathError <- rep(NA, nTestDays)
  for(k in 1:length(testDays)) {
    meanCaseError[k] <- mean(as.numeric(caseError[,k][caseError[,k] != Inf]), na.rm = T)
    meanDeathError[k] <- mean(as.numeric(deathError[,k][deathError[,k] != Inf]), na.rm = T)
  }

  deathErrNum <- caseErrNum <- matrix(NA, nrow = 50, ncol = nTestDays)
  for(k in 1:length(testDays)) {
    caseErrNum[,k] <- as.numeric(caseError[,k])
    deathErrNum[,k] <- as.numeric(deathError[,k])
  }
  caseErrNum[which(caseErrNum == Inf)] <- NA
  deathErrNum[which(deathErrNum == Inf)] <- NA
  
  return(list(caseError = caseError, deathError = deathError))
}


mase4_30 <- sirdPredError(as.Date("2020-05-01") + 0:20, "/Users/gregw/Dropbox/Projects/covidSird/Output/States/2020-08-30_4_30/Data/")

mase5_30 <- sirdPredError(as.Date("2020-05-31") + 0:20, "/Users/gregw/Dropbox/Projects/covidSird/Output/States/2020-08-30_5_30/Data/")

mase6_30 <- sirdPredError(as.Date("2020-07-01") + 0:20, "/Users/gregw/Dropbox/Projects/covidSird/Output/States/2020-08-30_6_30/Data/")

mase7_30 <- sirdPredError(as.Date("2020-07-31") + 0:20, "/Users/gregw/Dropbox/Projects/covidSird/Output/States/2020-08-30_7_30/Data/")


medianCaseError <- loCaseError <- upCaseError <-
  medianDeathError <- loDeathError <- upDeathError <-
  matrix(NA, nrow = 4, ncol = lagDays)

for(j in 1:lagDays) {
  medianCaseError[1,j] <- median(mase4_30$caseError[,j])
  loCaseError[1,j] <- quantile(mase4_30$caseError[,j], probs = .25)
  upCaseError[1,j] <- quantile(mase4_30$caseError[,j], probs = .75)
  medianCaseError[2,j] <- median(mase5_30$caseError[,j])
  loCaseError[2,j] <- quantile(mase5_30$caseError[,j], probs = .25)
  upCaseError[2,j] <- quantile(mase5_30$caseError[,j], probs = .75)
  medianCaseError[3,j] <- median(mase6_30$caseError[,j])
  loCaseError[3,j] <- quantile(mase6_30$caseError[,j], probs = .25)
  upCaseError[3,j] <- quantile(mase6_30$caseError[,j], probs = .75)
  medianCaseError[4,j] <- median(mase7_30$caseError[,j], na.rm = T)
  loCaseError[4,j] <- quantile(mase7_30$caseError[,j], probs = .25, na.rm = T)
  upCaseError[4,j] <- quantile(mase7_30$caseError[,j], probs = .75, na.rm = T)
  medianDeathError[1,j] <- median(mase4_30$deathError[,j])
  loDeathError[1,j] <- quantile(mase4_30$deathError[,j], probs = .25)
  upDeathError[1,j] <- quantile(mase4_30$deathError[,j], probs = .75)
  medianDeathError[2,j] <- median(mase5_30$deathError[,j])
  loDeathError[2,j] <- quantile(mase5_30$deathError[,j], probs = .25)
  upDeathError[2,j] <- quantile(mase5_30$deathError[,j], probs = .75)
  medianDeathError[3,j] <- median(mase6_30$deathError[,j])
  loDeathError[3,j] <- quantile(mase6_30$deathError[,j], probs = .25)
  upDeathError[3,j] <- quantile(mase6_30$deathError[,j], probs = .75)
  medianDeathError[4,j] <- median(mase7_30$deathError[,j], na.rm = T)
  loDeathError[4,j] <- quantile(mase7_30$deathError[,j], probs = .25, na.rm = T)
  upDeathError[4,j] <- quantile(mase7_30$deathError[,j], probs = .75, na.rm = T)
}

pdf(file = paste0("/Users/gregw/Dropbox/Projects/covidSird/Output/States/CaseMASE.pdf"), width = 9, height = 6)
par(mai = c(.8,.8,0,.2), mgp = c(3,.75,0))
plot(NA, xlim = c(1, lagDays), ylim = range(c(loCaseError, upCaseError), na.rm = T),
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n", bty = "n")
for(k in 1:4) {
  polygon(c(1:lagDays, lagDays:1), c(loCaseError[k,], rev(upCaseError[k,])), col = paste0(plotCols[k+1], "33"), border = F)
  points(1:lagDays, medianCaseError[k,], col = plotCols[k+1], type = "o", pch = 16)
}
axis(2, at = c(-9999999, 99999999), col = "grey33")
axis(1, at = c(-9999999, 99999999), col = "grey33")
axis(2, at = seq(0,5,.2), col = "grey33", las = 2, col.axis = "grey33",
labels = seq(0,5,.2), cex = .9, tick = F, hadj = .75)
axis(1, at = seq(7,lagDays,7), col.axis = "grey33", col = "grey33", tick = F, padj = -1.5)
 mtext("Mean Absolute Scaled Error", side=2, line=2.5, col="grey33", cex=1)
  mtext("Days from End of Training Data", side=1, line=1.5, col="grey33", cex=1)
dev.off()

pdf(file = paste0("/Users/gregw/Dropbox/Projects/covidSird/Output/States/DeathMASE.pdf"), width = 9, height = 6)
par(mai = c(.8,.8,0,.2), mgp = c(3,.75,0))
plot(NA, xlim = c(1, lagDays), ylim = range(c(loDeathError, upDeathError), na.rm = T),
     xlab = "", ylab = "",
     yaxt = "n", xaxt = "n", bty = "n")
for(k in 1:4) {
  polygon(c(1:lagDays, lagDays:1), c(loDeathError[k,], rev(upDeathError[k,])), col = paste0(plotCols[k+1], "33"), border = F)
  points(1:lagDays, medianDeathError[k,], col = plotCols[k+1], type = "o", pch = 16)
}
axis(2, at = c(-9999999, 99999999), col = "grey33")
axis(1, at = c(-9999999, 99999999), col = "grey33")
axis(2, at = seq(0,5,.2), col = "grey33", las = 2, col.axis = "grey33",
labels = seq(0,5,.2), cex = .9, tick = F, hadj = .75)
axis(1, at = seq(7,lagDays,7), col.axis = "grey33", col = "grey33", tick = F, padj = -1.5)
 mtext("Mean Absolute Scaled Error", side=2, line=2.5, col="grey33", cex=1)
  mtext("Days from End of Training Data", side=1, line=1.5, col="grey33", cex=1)
dev.off()





rbPal <- colorRampPalette(c("#5954A4", plotCols[c(3,2,1,4,6,5)]))
stCols <- rbPal(14)

# ---
pdf(file = paste0("/Users/gregw/Dropbox/COVID-19/Reports/2020-04-15 States/CaseMASE.pdf"), width = 9, height = 6)
par(mai = c(.8,.8,0,.2), mgp = c(3,.75,0))

boxplot(sqrt(caseError), outline = F, col = paste0(stCols, "40"), #paste0(c("#5954A4", plotCols[c(3,2,1,4,6,5)],"#5954A4", plotCols[c(3,2,1,4,6,5)]), c("40","80")), #paste0(c("#5954A4", plotCols[c(3,2,1,4,6,5)]), "40"),
        ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, 1.2), border = "grey33")
axis(2, at = seq(0,5,.2), col = "grey33", las = 2, col.axis = "grey33",
 labels = seq(0,5,.2), cex = .9, tick = F, hadj = .75)

 axis(1, at = 1:length(testDays), labels = gsub(" 0", " ", format(testDays, "%b %d")), col.axis = "grey33", col = "grey33", tick = F, padj = -2, cex.axis= .75)

abline(h = 1, col = "grey82", lwd = 2)
 mtext("Mean Absolute Scaled Error", side=2, line=2.5, col="grey33", cex=1)
dev.off()
# ---
pdf(file = paste0("/Users/gregw/Dropbox/COVID-19/Reports/2020-04-15 States/DeathMASE.pdf"), width = 9, height = 6)
par(mai = c(.8,.8,0,.2), mgp = c(3,.75,0))

boxplot(sqrt(deathError), outline = F, col = paste0(stCols, "40"),  #paste0(c("#5954A4", plotCols[c(3,2,1,4,6,5)]), "40"),
        ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, 1.2), border = "grey33")
axis(2, at = seq(0,5,.2), col = "grey33", las = 2, col.axis = "grey33",
 labels = seq(0,5,.2), cex = .9, tick = F, hadj = .75)

 axis(1, at = 1:length(testDays), labels = gsub(" 0", " ", format(testDays, "%b %d")), col.axis = "grey33", col = "grey33", tick = F, padj = -2, cex.axis= .75)
abline(h = 1, col = "grey82", lwd = 2)
 mtext("Mean Absolute Scaled Error", side=2, line=2.5, col="grey33", cex=1)
dev.off()
# ---




ny <- which(velocLogCases$loc == 1)
co <- which(velocLogCases$loc == 13)
wv <- which(velocLogCases$loc == 44)

options(scipen=999)
# raw cases

pdf(file = paste0("//Users/gregw/Dropbox/Projects/covidSird/Output/States/Cases.pdf"), width = 6, height = 6)
par(mai = c(.8,.8,0,.2))
plot((velocLogCases$u[ny]), col = plotCols[3], type = "l",
     ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n",
     lwd = 2, ylim = c(0, (max(velocLogCases$u[ny]))))
points((velocLogCases$u[co]), col = plotCols[2], type = "l", lwd = 2)
points((velocLogCases$u[wv]), col = plotCols[1], type = "l", lwd = 2)

axis(2, at = c(-9999999, 99999999), col = "grey33")
axis(1, at = c(-9999999, 99999999), col = "grey33")

axis(1, at = seq(0,200,50), col = "grey33", col.axis = "grey33", tick = F,
     padj = -1)
axis(2, at = seq(0,1000000,100000), labels = c("0", paste0(seq(100,1000,100), "k")), col = "grey33", col.axis = "grey33", tick = F, las = 2, hadj = .8)

 mtext("Cumulative Confirmed Cases", side=2, line=3, col="grey33", cex=1)
 mtext("Days since 100+ Confirmed Cases", side=1, line=2.5, col="grey33", cex=1)
dev.off()

# log cases

pdf(file = paste0("//Users/gregw/Dropbox/Projects/covidSird/Output/States/LogCases.pdf"), width = 6, height = 6)
par(mai = c(.8,.8,0,.2))
plot(log(velocLogCases$u[ny]), col = plotCols[3], type = "l",
     ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n",
     lwd = 2, ylim = c(0, log(max(velocLogCases$u[ny]))))
points(log(velocLogCases$u[co]), col = plotCols[2], type = "l", lwd = 2)
points(log(velocLogCases$u[wv]), col = plotCols[1], type = "l", lwd = 2)

axis(2, at = c(-9999999, 99999999), col = "grey33")
axis(1, at = c(-9999999, 99999999), col = "grey33")

axis(1, at = seq(0,200,50), col = "grey33", col.axis = "grey33", tick = F,
     padj = -1)
axis(2, at = log(c(1,10,100,1000,10000, 100000, 1000000)), labels = c("1","10","100","1,000","10,000", "100,000", "1,000,000"), col = "grey33", col.axis = "grey33", tick = F, las = 2, hadj = .8)

 mtext("Log Cumulative Confirmed Cases", side=2, line=3, col="grey33", cex=1)
 mtext("Days since 100+ Confirmed Cases", side=1, line=2.5, col="grey33", cex=1)
dev.off()

 # d/dt log cases
 pdf(file = paste0("//Users/gregw/Dropbox/Projects/covidSird/Output/States/VelocLogCases.pdf"), width = 6, height = 6)
 par(mai = c(.8,.8,0,.2))
 plot(lowess(velocLogCases$y[ny], f = .25), col = plotCols[3], type = "l",
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n",
      ylim = c(0, 0.25), lwd = 2)
 points(lowess(velocLogCases$y[co], f = .25), col = plotCols[2], type = "l", lwd = 2)
 points(lowess(velocLogCases$y[wv], f = .25), col = plotCols[1], type = "l", lwd = 2)

 axis(2, at = c(-9999999, 99999999), col = "grey33")
 axis(1, at = c(-9999999, 99999999), col = "grey33")

 axis(1, at = seq(0,200,50), col = "grey33", col.axis = "grey33", tick = F,
      padj = -1)
 axis(2, at = seq(0,.5,.05), col = "grey33", col.axis = "grey33", tick = F, las = 2, hadj = .5)

  mtext("Smoothed Velocity of Log Cumulative Confirmed Cases", side=2, line=2.5, col="grey33", cex=1)
  mtext("Days since 100+ Confirmed Cases", side=1, line=2.5, col="grey33", cex=1)
 dev.off()
