# NB: subject randomisation numbers have be over-written in the CSV file
#     to further anonymise the data

#**********************************************
# Read-in CSV data
#**********************************************

epileptic <- read.csv("./data-raw/qol/qol.csv")

colnames(epileptic)[1:4] <- c("id", "with.time", "trt", "with.status")

epileptic <- by(epileptic, epileptic$id, FUN = function(x) {
  fu <- x[rep(1, 4), 1:4]
  time <- as.numeric(x[1, c("time.b", "time.3m", "time.1y", "time.2y")])
  anxiety <- as.numeric(x[1, c("Anxiety.b", "Anxiety.3m", "Anxiety.1y", "Anxiety.2y")])
  depress <- as.numeric(x[1, c("Depress.b", "Depress.3m", "Depress.1y", "Depress.2y")])
  aep <- as.numeric(x[1, c("AEP.b", "AEP.3m", "AEP.1y", "AEP.2y")])
  cbind(fu, time, anxiety, depress, aep)
})

#------------------------------------------------------------------

#**********************************************
# Code missing data
#**********************************************

epileptic <- do.call("rbind", epileptic)
epileptic$time[epileptic$time == "999"] <- NA
epileptic$time[epileptic$time < 0] <- NA
epileptic$anxiety[epileptic$anxiety == "999"] <- NA
epileptic$depress[epileptic$depress == "999"] <- NA
epileptic$aep[epileptic$aep == "999"] <- NA

epileptic <- epileptic[!is.na(epileptic$time), ]

#------------------------------------------------------------------

#**********************************************
# Finalise datasets
#**********************************************

epileptic$with.status2 <- as.numeric(epileptic$with.status != 0)
row.names(epileptic) <- seq_len(nrow(epileptic))
epileptic.qol <- epileptic
epileptic.qol <- droplevels(qol)
save(epileptic.qol, file = "epileptic.qol.Rdata")
