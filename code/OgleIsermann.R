################################################################################
################################################################################
##
## Analysis script for ... Ogle, DH.  201X.  .  Fisheries
##   Research XX:XXX-XXX.
##
## Need to be patient with bootstrapping functions.
## Expected that data files are in a "data" directory within ghe current
##   working directory.
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
library(plotrix)

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
##
##   1. Fit all models
##   2. Estimate t480 by rearranging typical and original VBGF
##   3. Estimate t480 as parameter in new VBGF
##   4. Make summary table to show equivalencies
##   5. Assessed alternative starting values
##
################################################################################

## Load data ... and isolate BBN
df <- read.csv("data/LMWhitefish_byStock.csv")
bbn <- droplevels(subset(df,Stock=="BBN"))

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
                    predict(bootT,vbT,t=ages,conf.level=cl)[,3:4]),
              RSS=c(sum(residuals(fitT)^2),NA,NA))
rownames(resT)[4:6] <- c("tr",paste0("pred",ages))

resO <- rbind(cbind(Ests=coef(fitO),confint(bootO,conf.level=cl)),
              c(trO,predict(bootO,calc_trO,Lr=Lr,conf.level=cl)[,3:4]),
              cbind(vbO(ages,coef(fitO)),
                    predict(bootO,vbO,t=ages,conf.level=cl)[,3:4]),
              RSS=c(sum(residuals(fitO)^2),NA,NA))
rownames(resO)[4:6] <- c("tr",paste0("pred",ages))

bootN2 <- bootN
bootN2$coefboot <- cbind(bootN2$coefboot,Lr)
resN <- rbind(cbind(Ests=coef(fitN,Lr),confint(bootN,conf.level=cl)),
              cbind(vbnew(ages,c(coef(fitN),Lr)),
                    predict(bootN2,vbnew,t=ages,conf.level=cl)[,3:4]),
              RSS=c(sum(residuals(fitN)^2),NA,NA))
rownames(resN)[4:5] <- c(paste0("pred",ages))

round(resT,4)
round(resO,4)
round(resN,4)
aictab(list(fitT,fitO,fitN),c("Typical","Original","New"))

# Alternative starting values .. Just checking for convergence and relationship
#   to parameter estimates from the fits above.
svT1 <- list(Linf=700,K=0.2,t0=-2)
fitT1 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT1)
svT2 <- list(Linf=400,K=0.5,t0=2)
fitT2 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT2)
svT3 <- list(Linf=650,K=0.1,t0=-2)
fitT3 <- nls(Len~vbT(Age,Linf,K,t0),data=bbn,start=svT3)
round(cbind(coef(fitT),coef(fitT1),coef(fitT2),coef(fitT3)),3)    # OK

svO1 <- list(Linf=700,K=0.2,L0=100)
fitO1 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO1)
svO2 <- list(Linf=400,K=0.2,L0=300)
fitO2 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO2)
svO3 <- list(Linf=650,K=0.1,L0=100)
fitO3 <- nls(Len~vbO(Age,Linf,K,L0),data=bbn,start=svO3)
round(cbind(coef(fitO),coef(fitO1),coef(fitO2),coef(fitO3)),3)    # OK

svN1 <- list(Linf=700,K=0.5,tr=6)
fitN1 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN1)
svN2 <- list(Linf=400,K=0.2,tr=10)
fitN2 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN2)
svN3 <- list(Linf=650,K=0.1,tr=10)
fitN3 <- nls(Len~vbnew(Age,Linf,K,tr,Lr),data=bbn,start=svN3)
round(cbind(coef(fitN),coef(fitN1),coef(fitN2),coef(fitN3)),3)    # OK






################################################################################
################################################################################
## Lake Winnibigoshis Walleye analysis to compare sexes within a year
##
##   1. Fit all models
##   2. Assess support for each model
##   3. Make summary table
##   4. Make summary graphics
##
################################################################################

# Load data
df <- read.csv("data/WBGWalleye.csv",na.strings=".")
wae12 <- droplevels(subset(df,Year==2012 & !is.na(Sex)))

# Set critical length
Lr <- round(17*25.4,0)

# Create all models
vbOm   <- TL~Lr+(Linf-Lr)*(1-exp(-K*(Age-tr)))
vbLKTr <- TL~Lr+(Linf[Sex]-Lr)*(1-exp(-K[Sex]*(Age-tr[Sex])))
vbLK   <- TL~Lr+(Linf[Sex]-Lr)*(1-exp(-K[Sex]*(Age-tr)))
vbLTr  <- TL~Lr+(Linf[Sex]-Lr)*(1-exp(-K*(Age-tr[Sex])))
vbKTr  <- TL~Lr+(Linf-Lr)*(1-exp(-K[Sex]*(Age-tr[Sex])))
vbL    <- TL~Lr+(Linf[Sex]-Lr)*(1-exp(-K*(Age-tr)))
vbK    <- TL~Lr+(Linf-Lr)*(1-exp(-K[Sex]*(Age-tr)))
vbTr   <- TL~Lr+(Linf-Lr)*(1-exp(-K*(Age-tr[Sex])))

# Create all starting values
svOm   <- list(Linf=650,K=0.2,tr=2)
svLKTr <- Map(rep,svOm,c(2,2,2))
svLK   <- Map(rep,svOm,c(2,2,1))
svLTr  <- Map(rep,svOm,c(2,1,2))
svKTr  <- Map(rep,svOm,c(1,2,2))
svL    <- Map(rep,svOm,c(2,1,1))
svK    <- Map(rep,svOm,c(1,2,1))
svTr   <- Map(rep,svOm,c(1,1,2))

# Fit all models
fitOm   <- nls(vbOm,  data=wae12,start=svOm)
fitLKTr <- nls(vbLKTr,data=wae12,start=svLKTr)
fitLK   <- nls(vbLK,  data=wae12,start=svLK)
fitLTr  <- nls(vbLTr, data=wae12,start=svLTr)
fitKTr  <- nls(vbKTr, data=wae12,start=svKTr)
fitL    <- nls(vbL,   data=wae12,start=svL)
fitK    <- nls(vbK,   data=wae12,start=svK)
fitTr   <- nls(vbTr,  data=wae12,start=svTr)

extraSS(fitOm,com=fitLKTr,sim.name="{Omega}",com.name="{L,K,tr}")
extraSS(fitLK,fitLTr,fitKTr,sim.name=c("{Linf,K}","{Linf,tr}","{K,tr}"),
        com=fitLKTr,com.name="{Linf,K,tr}")
extraSS(fitL,fitTr,sim.name=c("{Linf}","{tr}"),com=fitLTr,com.name="{Linf,tr}")

summary(fitLTr,correlation=TRUE)
( resLTr <- rbind(cbind(Ests=coef(fitLTr),confint(fitLTr))) )

## Figure 2
pdf("results/Figure_2.PDF",width=4,height=4)
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2.1,0.4,0),tcl=-0.2,las=1,xaxs="i",yaxs="i")
jit <- 0.1; pt.cex <- 0.8; clr <- col2rgbt("black",1/2)
ltys <- 3:2; pchs <- 0:1; sfrac <- 0.005; slwd=0.7
# base plot
plot(TL~Age,data=wae12,xlab="Age (years)",ylab="Total Length (mm)",
     xlim=c(0,13.1),ylim=c(200,710),col="white")
axis(1,0:13,NA,tcl=-0.1)
# show Lr
axis(2,Lr,expression(italic(L)[r]))
abline(h=Lr,col="gray90")
# show points
points(TL~I(Age-jit),data=subset(wae12,Sex=="F"),pch=pchs[1],col=clr,cex=pt.cex^2)
points(TL~I(Age+jit),data=subset(wae12,Sex=="M"),pch=pchs[2],col=clr,cex=pt.cex)
# show best-fit lines
curve(vbnew(x,c(coef(fitLTr)[c(1,3,4)],Lr)),from=0,to=18,add=TRUE,
      lwd=2,lty=ltys[1])
curve(vbnew(x,c(coef(fitLTr)[c(2,3,5)],Lr)),from=0,to=18,add=TRUE,
      lwd=2,lty=ltys[2])
# Show Linf
abline(h=resLTr[c("Linf1","Linf2"),"Ests"],lty=ltys,lwd=1)
plotCI(rep(0.25,2),resLTr[c("Linf1","Linf2"),"Ests"],
       li=resLTr[c("Linf1","Linf2"),"2.5%"],
       ui=resLTr[c("Linf1","Linf2"),"97.5%"],add=TRUE,pch=NA,
       sfrac=sfrac,lwd=slwd)
# Show tr
lines(rep(resLTr["tr1","Ests"],2),c(Lr,0),lty=ltys[1],lwd=1)
lines(rep(resLTr["tr2","Ests"],2),c(Lr,0),lty=ltys[2],lwd=1)
plotCI(resLTr[c("tr1","tr2"),"Ests"],rep(210,2),
       li=resLTr[c("tr1","tr2"),"2.5%"],
       ui=resLTr[c("tr1","tr2"),"97.5%"],err="x",add=TRUE,pch=NA,
       sfrac=sfrac,lwd=slwd)
dev.off()
