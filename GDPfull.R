#library(seasonal)
## Data for the monthly truth and quarterly response 
fulldata <- read.csv("C:/Users/mosleyl/OneDrive - Lancaster University/PhD/SparseTD/Paper 1 - GDP data and sort out/BigdataGDP.csv",header = T)
quarterlyGDP <- read.csv("C:/Users/mosleyl/OneDrive - Lancaster University/PhD/SparseTD/Paper 1 - GDP data and sort out/quarterlyCVMSAdata.csv",header=T)

GDPmonthlyinx <- fulldata$GDPM[1:150]
GDPmonthlyinx_TS <- ts(GDPmonthlyinx, deltat = 1/12, start = c(2008,1))
GDPquarterly_CVMSA_TS <- ts(quarterlyGDP$CVMSA[1:50], deltat = 1/4, start = c(2008,1))

m2 <- td(GDPquarterly_CVMSA_TS ~ 0 + GDPmonthlyinx_TS, method = "denton-cholette")
GDPmonthlyPRICE <- m2$values

scaleY <- scale(GDPquarterly_CVMSA_TS)
quarterscale_center <- attr(scaleY,"scaled:center")
quarterscale_scale <- attr(scaleY,"scaled:scale")
tsY <- ts(scaleY, deltat = 1/4, start = c(2008,1))

## Indicator series
# Aggregate
aggdata <- read.csv("C:/Users/mosleyl/OneDrive - Lancaster University/PhD/SparseTD/Paper 1 - GDP data and sort out/Project1Data.csv",header = T)
newaggregates <- read.csv("C:/Users/mosleyl/OneDrive - Lancaster University/PhD/SparseTD/Paper 1 - GDP data and sort out/newaggregates.csv",header = T)
aggdata2 <- aggdata[,-c(1,10,11)]
aggdata2[,c(17,18)] <- newaggregates
X_NSA <- data.matrix(aggdata2)
X_aggagg_NSA <- X_NSA[,c(1:8,17,18)]

matrixSA <- function(mat) {
  SAX <- matrix(NA, ncol = ncol(mat), nrow(mat))
  for(i in 1:ncol(mat)) {
    ind <- ts(mat[,i], deltat = 1/12, start = c(2008,1))
    SAind <- seas(ind)
    SA <- as.numeric(SAind$data[,1])
    print("done")
    SAX[,i] <- SA
  }
  return(SAX)
}

X_aggagg_SA <- matrixSA(X_aggagg_NSA)

tsXagg <- ts(scale(X_aggagg_SA), deltat = 1/12, start = c(2008,1))

# Full 97
fulldata <- fulldata[,-c(76,81)]
fulldata$NSA <- na.interpolation(fulldata$NSA)
fulldata$NSA.1 <- na.interpolation(fulldata$NSA.1)
fulldata$NSA.2 <- na.interpolation(fulldata$NSA.2)
fulldata$NSA.3 <- na.interpolation(fulldata$NSA.3)
fulldata$NSA.5 <- na.interpolation(fulldata$NSA.5)
fulldata$NSA.6 <- na.interpolation(fulldata$NSA.6)
fulldata$NSA.7 <- na.interpolation(fulldata$NSA.7)
fulldata$NSA.8 <- na.interpolation(fulldata$NSA.8)
fulldata <- fulldata[,-1]
fulldata <- fulldata[,-c(39,40,41)]
Xfull_NSA <- data.matrix(fulldata)
xfull_SA2 <- matrixSA(Xfull_NSA)
#Xfull_SA3 <- matrixSA(Xfull_NSA) # see below
tsXfull <- ts(scale(Xfull_SA4), deltat = 1/12, start = c(2008,1))

# CL applied
m3 <- td(tsY ~ 0 + tsXagg, to = "month", truncated.rho = 0)
CLest <- m3$values
unscaledCL2 <- CLest*quarterscale_scale + quarterscale_center/3    # using est*SD(quarterly)+mean(quarterly)/3

# spTD applied 
m4 <- SparseTD.estimates(tsXagg, tsY)
SPTDest <- ts(m4$rfseries.estimate, deltat = 1/12, start = c(2008,1))
unscaledSPTD <- SPTDest*quarterscale_scale + quarterscale_center/3 # unscaling the estimates 

# spTD 97 applied 
m5 <- SparseTD.estimates(tsXfull,tsY)
SPTD97est <- ts(m5$rfseries.estimate, deltat = 1/12, start = c(2008,1))
unscaledSPTD97 <- SPTD97est*quarterscale_scale + quarterscale_center/3 # unscaling the estimates 






