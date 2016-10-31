################################################################################
################################################################################
##
## Analysis script for ... Ogle, DH.  201X.  .  Fisheries
##   Research XX:XXX-XXX.
##
## Need to be patient with bootstrapping functions.
## May need to create a directory called "results" in your current working
##   directory to hold the figures produced by pdf() (or not run the pdf()
##   and dev.off() functions to simply produce the figures on the local device).
##   Could use (in R) to create the directory (assumes that you have set your
##   working directory to the same location as this script) ...
##
##   dir.create("results")
##
## This code was tested on a Windows 7 machine using 32-bit R v3.3.1 and a
##   Macintosh (El Capitan OS) machine using 64-bit R v3.3.1.  The code runs
##   without error on both machines.
##
################################################################################
################################################################################

################################################################################
## SETUP
## 
## Requires FSA (>=0.8.11) from CRAN, installed with:
##
##   install.packages("FSA")
##
################################################################################
## Load required packages
library(FSA)
library(nlstools)
library(AICcmodavg)

## Set random number seed so that bootstrappings are repeatable
set.seed(347834)

## Create a function with the new, typical, and original VBGFs
vbnew <- vbFuns("Ogle")
vbT <- vbFuns("Typical")
vbO <- vbFuns("Original")

## Create a function to compute tr from typical VBGF results
calc_trT <- function(Lr,Linf,K=NULL,t0=NULL) {
  if (length(Linf)==3) {
    K <- Linf[[2]]
    t0 <- Linf[[3]]
    Linf <- Linf[[1]] }
  (log(1-Lr/Linf))/(-K)+t0
}

## Create a function to compute tr from original VBGF results
calc_trO <- function(Lr,Linf,K=NULL,L0=NULL) {
  if (length(Linf)==3) {
    K <- Linf[[2]]
    L0 <- Linf[[3]]
    Linf <- Linf[[1]] }
  (log((Linf-Lr)/(Linf-L0)))/(-K)
}

## Load data ... and isolate BBN
df <- read.csv("data/LMWhitefish2.csv")
str(df)
bbn <- droplevels(subset(df,Zone=="WI-2"))

################################################################################
################################################################################
## Figure 1
################################################################################
pdf("results/Figure_1.PDF",width=4,height=4)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2.1,0.4,0),tcl=-0.2,las=1,xaxs="i",yaxs="i")
Linf <- 250; t0 <- -0.7; K <- 0.5
curve(vbnew(x,Linf,K,t0,0),from=-1,to=5,lwd=3,
      ylab="Length",xlab="Age",ylim=c(0,1.05*Linf),xlim=c(-1,5),yaxt="n",xaxt="n")
axis(1,-1:5,c(NA,0:5))
axis(2,seq(0,250,50),c(0,NA,100,NA,200,NA))
# Mark Linf on plot
abline(h=Linf,lwd=2,lty=2,col="gray50")
axis(2,Linf,expression(italic(L)[infinity]),cex.axis=1.25,lwd=2,tick=FALSE)
# Mark t0 on plot
axis(1,t0,expression(italic(t)[0]),cex.axis=1.25,padj=0.3,tick=FALSE)
points(t0,0,cex=1.2,xpd=TRUE,pch=21,col="black",bg="gray80")
# Mark L0 on plot
L0 <- vbnew(0,Linf,K,t0,0)
axis(2,L0,expression(italic(L)[0]),cex.axis=1.25,padj=0.3,tick=FALSE)
lines(c(0,0),c(0,L0),lwd=2,lty=3,col="gray50")
lines(c(-2,0),c(L0,L0),lwd=2,lty=3,col="gray50")
points(0,L0,cex=1.2,xpd=TRUE,pch=21,col="black",bg="gray80")
# Mark an (tr,Lr) on plot
tr <- 1.5
Lr <- vbnew(tr,Linf,K,t0,0)
axis(2,Lr,expression(italic(L)[r]),cex.axis=1.25,padj=0.3,tick=FALSE)
axis(1,tr,expression(italic(t)[r]),cex.axis=1.25,padj=0.3,tick=FALSE)
lines(c(tr,tr),c(0,Lr),lwd=2,lty=3,col="gray50")
lines(c(-2,tr),c(Lr,Lr),lwd=2,lty=3,col="gray50")
points(tr,Lr,cex=1.2,xpd=TRUE,pch=21,col="black",bg="gray80")
dev.off()


################################################################################
################################################################################
## Big Bay de Noc analysis
##   1. Fit all models
##   2. Estimate t480 by rearranging typical and original VBGF
##   3. Estimate t480 as parameter in new VBGF
##   4. Make summary table to show equivalencies
##
##   Note that alternative starting values are further below.
################################################################################
# Set critical length
Lr <- 480
# Set ages at which to make some predictions
ages <- c(8,20)

# Fit the three models
svT <- list(Linf=550,K=0.3,t0=0)
fitT <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT)
bootT <- nlsBoot(fitT)
svO <- list(Linf=550,K=0.3,L0=200)
fitO <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO)
bootO <- nlsBoot(fitO)
svN <- list(Linf=550,K=0.3,tr=8)
fitN <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN)
bootN <- nlsBoot(fitN)

# Compute tr estimates from each model
trT <- calc_trT(Lr,coef(fitT))
trO <- calc_trO(Lr,coef(fitO))
trN <- coef(fitN)[["tr"]]

# Put some results together to make Table 1
cl <- 0.90  # confidence level
resT <- rbind(cbind(Ests=coef(fitT),confint(bootT,conf.level=cl)),
              c(trT,predict(bootT,calc_trT,Lr=Lr,conf.level=cl)[,3:4]),
              cbind(vbT(ages,coef(fitT)),
                    predict(bootT,vbT,t=ages,conf.level=cl)[,3:4]))
rownames(resT)[4:6] <- c("tr",paste0("pred",ages))

resO <- rbind(cbind(Ests=coef(fitO),confint(bootO,conf.level=cl)),
              c(trO,predict(bootO,calc_trO,Lr=Lr,conf.level=cl)[,3:4]),
              cbind(vbO(ages,coef(fitO)),
                    predict(bootO,vbO,t=ages,conf.level=cl)[,3:4]))
rownames(resO)[4:6] <- c("tr",paste0("pred",ages))

bootN2 <- bootN
bootN2$coefboot <- cbind(bootN2$coefboot,480)
resN <- rbind(cbind(Ests=coef(fitN,Lr),confint(bootN,conf.level=cl)),
              cbind(vbnew(ages,c(coef(fitN),Lr)),
                    predict(bootN2,vbnew,t=ages,conf.level=cl)[,3:4]))
rownames(resN)[4:5] <- c(paste0("pred",ages))

round(resT,4)
round(resO,4)
round(resN,4)
aictab(list(fitT,fitO,fitN),c("Typical","Original","New"))



################################################################################
################################################################################
## Big Bay de Noc and Green Bay analysis
##   1. Fit all models
##   2. Estimate t480 by rearranging typical and original VBGF
##   3. Estimate t480 as parameter in new VBGF
##   4. Make summary table to show equivalencies
##
##   Note that alternative starting values are further below.
################################################################################
# Set critical length
Lr <- 480
# Create all models
vbOm   <- Len~Lr+(Linf-Lr)*(1-exp(-K*(Age-tr)))
vbLKTr <- Len~Lr+(Linf[Zone]-Lr)*(1-exp(-K[Zone]*(Age-tr[Zone])))
vbLK   <- Len~Lr+(Linf[Zone]-Lr)*(1-exp(-K[Zone]*(Age-tr)))
vbLTr  <- Len~Lr+(Linf[Zone]-Lr)*(1-exp(-K*(Age-tr[Zone])))
vbKTr  <- Len~Lr+(Linf-Lr)*(1-exp(-K[Zone]*(Age-tr[Zone])))
vbL    <- Len~Lr+(Linf[Zone]-Lr)*(1-exp(-K*(Age-tr)))
vbK    <- Len~Lr+(Linf-Lr)*(1-exp(-K[Zone]*(Age-tr)))
vbTr   <- Len~Lr+(Linf-Lr)*(1-exp(-K*(Age-tr[Zone])))
# Create all starting values
svOm   <- list(Linf=550,K=0.25,tr=10)
svLKTr <- Map(rep,svOm,c(2,2,2))
svLK   <- Map(rep,svOm,c(2,2,1))
svLTr  <- Map(rep,svOm,c(2,1,2))
svKTr  <- Map(rep,svOm,c(1,2,2))
svL    <- Map(rep,svOm,c(2,1,1))
svK    <- Map(rep,svOm,c(1,2,1))
svTr   <- Map(rep,svOm,c(1,1,2))
# Fit all models
fitOm   <- nls(vbOm,  data=df,start=svOm)
fitLKTr <- nls(vbLKTr,data=df,start=svLKTr)
fitLK   <- nls(vbLK,  data=df,start=svLK)
fitLTr  <- nls(vbLTr, data=df,start=svLTr)
fitKTr  <- nls(vbKTr, data=df,start=svKTr)
fitL    <- nls(vbL,   data=df,start=svL)
fitK    <- nls(vbK,   data=df,start=svK)
fitTr   <- nls(vbTr,  data=df,start=svTr)
## Put together AIC results
aictab(list(fitLKTr,fitLK,fitLTr,fitKTr,fitL,fitK,fitTr,fitOm),
       c("LKTr","LK","LTr","KTr","L","K","Tr","Omega"))




################################################################################
################################################################################
## Testing different starting values for model fits with the Big Bay de Noc data
##   Just checking for convergence and relationship to parameter estimates
##   from the fits above.
################################################################################
################################################################################
svT1 <- list(Linf=700,K=0.2,t0=-2)
fitT1 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT1)
svT2 <- list(Linf=400,K=0.5,t0=2)
fitT2 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT2)
svT3 <- list(Linf=550,K=0.1,t0=-2)
fitT3 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT3)
round(cbind(coef(fitT),coef(fitT1),coef(fitT2),coef(fitT3)),3)    # OK

svO1 <- list(Linf=700,K=0.2,L0=100)
fitO1 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO1)
svO2 <- list(Linf=400,K=0.3,L0=300)
fitO2 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO2)
svO3 <- list(Linf=550,K=0.1,L0=100)
fitO3 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO3)
round(cbind(coef(fitO),coef(fitO1),coef(fitO2),coef(fitO3)),3)    # OK

svN1 <- list(Linf=700,K=0.5,tr=6)
fitN1 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN1)
svN2 <- list(Linf=400,K=0.3,tr=10)
fitN2 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN2)
svN3 <- list(Linf=550,K=0.1,tr=10)
fitN3 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN3)
round(cbind(coef(fitN),coef(fitN1),coef(fitN2),coef(fitN3)),3)    # OK
