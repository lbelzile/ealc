#Data analysis of ``Extremal attractors of Liouville copulas''
#(C) Leo Belzile
#Most routines are found in the packages ExtLiouv, with some utilities from
#packages "mev" and "lcopula"
#All of the packages are available from the Github page of the first author
rm(list=ls())
graphics.off()
#Add path to user directory
pathtofile <- ""
setwd(pathtofile)
if(!"Danube.RData" %in% list.files()){
 warning("User must provide a valid working directory containing `Danube.RData'")
}
devtools::install_github("lbelzile/mev")
#Extreme value and time series libraries, can be installed from CRAN
#install.packages("evd", dependencies=TRUE)
library(evd);  library(ismev); library(xts); library(lubridate)
#Libraries used in the analysis, whose up-to-date versions are available viz
#devtools::install_github("lbelzile/mev")
library(ExtLiouv); library(mev);
#Library to export images and further options
library(tikzDevice)
options(tikzLatexPackages =
          c("\\usepackage{tikz}\n",
            "\\usepackage[active,tightpage,psfixbb]{preview}\n",
            "\\PreviewEnvironment{pgfpicture}\n",
            "\\setlength\\PreviewBorder{0pt}\n",
            "\\usepackage{times}\n"
          ))
setTikzDefaults(overwrite = FALSE)

###############################################################
####          Data analysis based on Danube data           ####
###############################################################

#Load Danube data
load("Danube.RData")
plots <- TRUE
### Preliminary analysis
#Extract data, format to time series
#Extract three stations on the Isar river
wst <- c(16,17,19)
d <- length(wst)
dat <- zoo::na.trim(river[,wst],"both", is.na="any")
#Any missing values: no
summary(dat)

#How many years of data?
time(head(dat[,1],1)); time(tail(dat[,1],1))
colnames(dat) <- c("Muenchen","Pupplinger Au","Lenggries")
if(plots){
  #Plot the time series
  tikz("figure/Isar_river_flow.tex", standAlone=TRUE, width=8, height=5)
  plot.xts(dat, multi.panel=TRUE,major.ticks = "year", minor.ticks = NULL, col=c(1,1,1), grid.col="lightgray",
           grid.ticks.on = "year", yaxis.same = FALSE, format.labels = "%Y",main="Isar daily river flow")
  dev.off()
}
#Concordance
cor(dat, method="spearman")
#Retain only summer floods, June to August
dat <- dat[month(dat) %in% 6:8,]

#Extract station coordinates
stations <-  coord[wst,]
N <- nrow(dat)

if(plots){
  #Threshold selection
  for(i in 1:d){
    mev::W.diag(dat[,i], q1 = 0.8, q2 = 0.99, model="nhpp", k=19, plots="PS", UseQuantiles = FALSE,)
    mev::NC.diag(x=as.vector(unlist(dat[,i])), GP.fit="ismev", u=quantile(dat[,i], seq(0.8,0.98,by=0.02)), size=0.05)
  }
}
#Settle for
qu <- rep(0.92,d)

if(plots){
  #Extremal index and serial dependence
  par(mfrow=c(2,2))
  for(i in 1:d){
    mev::ext.index(dat[,i], seq(0.9,0.997,by=0.001), c("wls","mle","inter"),TRUE)
  }
}

#Index of exceedances
exceedances.index <- isAbove(dat, threshold = (u <- sapply(1:d, function(i){
  quantile(dat[,i], probs=qu[i])})))
univ.exceed.index <- which(rowSums(exceedances.index)>=1)
#Create a matrix for the cluster maxima, shrinking later if necessary
exceedances <- matrix(nrow=length(univ.exceed.index),ncol=d)
for(i in 1:nrow(exceedances)){
  se <- max((univ.exceed.index[i]-3),1):min((univ.exceed.index[i]+3))
  for(j in 1:d){
    exceedances[i,j] <- max(dat[se,j])
  }
  if(i>1){
    if(any(apply(exceedances[(i-1):i,],2, duplicated))){
      exceedances[(i-1):i,] <- t(matrix(rep(apply(exceedances[(i-1):i,], 2, max),2),ncol=2))
    }
  }
}
#Keep the non-duplicate
x <- unique(exceedances)
x[which(rowSums(apply(x,2, duplicated))>=1),]
#These duplicates are due after analysis to discreteness of measurements
n <- nrow(x)

#Univariate peaks-over-thresholds
#Need to deal with discreteness of the observations
#precision of 1, 0.1 and 0,1
fit.pgpd <- function(par, precis, vobs, u){
  vobs <- vobs[(vobs-precis/2-u)>0]
  -sum(log(pgpd(q = vobs+precis/2, loc = u, scale = exp(par[1]), shape = par[2])-
             pgpd(q = vobs-precis/2, loc = u, scale = exp(par[1]), shape = par[2])))
}

#Fit the Generalized Pareto
param.gpd <- function(i){
  fit <- mev::gp.fit(x[,i],u[i])$est
  fit2 <- optim(par=c(log(fit[1]),0.01), fn=fit.pgpd, u=u[i],
                vobs=x[,i],precis=c(1,0.1,0.1)[i], method="Nelder-Mead",control = list(reltol=1e-10))
  return(c(exp(fit2$par[1]),fit2$par[2]))
}
marg_par <- sapply(1:d,param.gpd)
marg_par_cont <- sapply(1:d, function(i){mev::gp.fit(x[,i],u[i])$est})
#No discernible difference, impact of rounding is negligeable
scale <- marg_par[1,]; shape <- marg_par[2,]
names(scale) <- sapply(1:d, function(i){ paste0("scale",i)})
names(shape) <- sapply(1:d, function(i){ paste0("shape",i)})
scale; shape


#Fit extreme value model - pairwise censored likelihood of Ledford and Tawn
#Scaled extremal Dirichlet model
fit1a <- fmvcpot(dat=x, u=u,  rho=1, model="dir",
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)

fit3a <- fmvcpot(dat=x, u=u,  rho=1,  model="dir", scale=scale, shape=shape, rho=1,
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit4a <- fmvcpot(dat=x, u=u, model="dir", scale=scale, shape=shape,
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit2a <- fmvcpot(dat=x, u=u, model="dir",
                 start=list(scale=log(marg_par[1,]), shape=marg_par[2,],
                            alpha=log(fit4a$estimate$alpha),rho=log(fit4a$estimate$rho)),
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="BFGS", transform=TRUE)
fit5a <- fmvcpot(dat=x, u=u,  alpha=c(1,1,1), model="dir",
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="BFGS", transform=TRUE)
fit6a <- fmvcpot(dat=x, u=u,  model="dir",
                 cscale=FALSE, cshape=TRUE, std.err=TRUE, method="BFGS", transform=TRUE)
#This fit diverges to the boundary of the parameter space, and does not converge with Nelder-Mead

#Scaled extremal Dirichlet model
fit1b <- fmvcpot(dat=x, u=u,  rho=1, model="negdir",
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit3b <- fmvcpot(dat=x, u=u,  rho=1,  model="negdir", scale=marg_par[1,], shape=marg_par[2,], rho=1,
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit4b <- fmvcpot(dat=x, u=u, model="negdir", scale=marg_par[1,], shape=marg_par[2,],
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit2b <- fmvcpot(dat=x, u=u, model="negdir", cscale=FALSE, cshape=FALSE,
                 start=list(scale=log(marg_par[1,]), shape=marg_par[2,],
                            alpha=log(fit4b$estimate$alpha),rho=log(0.35)),
                 std.err=TRUE, method="Nelder-Mead", transform=TRUE)
fit5b <- fmvcpot(dat=x, u=u,  alpha=c(1,1,1), model="negdir",
                 cscale=FALSE, cshape=FALSE, std.err=TRUE, method="BFGS", transform=TRUE)
fit6b <- fmvcpot(dat=x, u=u,  model="negdir",
                 cscale=FALSE, cshape=TRUE, std.err=TRUE, method="BFGS", transform=TRUE)

##Different fits

fit1a$par; fit1b$par #1: with rho=1 (Coles-Tawn model)
fit2a$par; fit2b$par #2: full model, no restriction
fit3a$par; fit3b$par #3: two stage procedure, rho=1
fit4a$par; fit4b$par #4: two stage procedure
fit5a$par; fit5b$par #5: (negative) logistic
fit6a$par; fit6b$par #6: same shape parameters

#This step is computationally intensive
#b <- compositemat(x, fitted=fit2b, B=1000, use.start=TRUE)
#b <- list(sensitivity=Hmat, godambe=godambe, variability=Jmat, parreplicates=bootpar)
#save(b = b, file="code/boot_composite.RData")
load("code/boot_composite.RData")
Hmat <- b$sensitivity

#Composite likelihood ratio test
#Test whether the reduction to the Dirichlet model
CLRT.stat <- 2*(fit1a$value-fit2b$value)
weight.CLRT <- solve(b$godambe)[10,10]/solve(Hmat)[10,10]
1-pchisq(CLRT.stat/weight.CLRT,1)
#Difference is significant

#Test whether we can reduce the model to the negative logistic model
CLRT.stat2 <- 2*(fit5a$value-fit2b$value)
weight.CLRT <- eigen(solve(solve(Hmat)[7:9,7:9])%*%solve(b$godambe)[7:9,7:9],only.values = TRUE)$values
1-momentchi2::sw(x=CLRT.stat2, weight.CLRT)
1-momentchi2::hbe(x=CLRT.stat2, weight.CLRT)
#Difference not significant

#Test whether we can restrict the model to the logistic
CLRT.stat3 <- 2*(fit5b$value-fit2b$value)
1-momentchi2::sw(x=CLRT.stat3, weight.CLRT)
1-momentchi2::hbe(x=CLRT.stat3, weight.CLRT)
#Difference not significant

#Quantile-quantile plot for the generalized Pareto distribution
#With pointwise 0.95 confidence interval based on Beta distribution of order statistics
#Adapted from evd:::qq.pot

gp.qqplot <- function (x, gp.obj, main = "Quantile-quantile plot", xlab = "Theoretical quantiles",
                       ylab = "Sample quantiles", ...){
  if(is.matrix(x)){
    stop("Data provided must be a vector, not a matrix")
  }
  dat <- sort(x[x>gp.obj$threshold])
  n <- length(dat)
  confint_lim <- t(sapply(1:n, function(i){
    qgpd(qbeta(c(0.025,0.975),i, n-i+1),loc=gp.obj$threshold,
         scale=gp.obj[['estimate']][1], shape=gp.obj[['estimate']][2])}))
  quant <- qgpd(rank(dat)/(n+1), loc = gp.obj$threshold, scale = gp.obj[['estimate']][1],
                shape = gp.obj[['estimate']][2])
  matplot(quant, cbind(dat, confint_lim), main = main, xlab = xlab,
          ylab = ylab, type = "pll", pch = 20,col=c(1,"grey","grey"),lty=c(1,2,2),bty="l",pty="s",...)
  abline(0, 1)
  invisible(cbind(quant,confint_lim))
}


#Marginal QQ-plots
graphics.off()
par(mfrow=c(1,3))
#Marginal plots, with estimates
for(i in 1:d){
  gp.qqplot(x=x[,i], gp.obj = list(threshold=u[i], estimate=c(scale[i], shape[i])))
}
#Marginal estimates for scaled Dirichlet
for(i in 1:d){
  gp.qqplot(x=x[,i], gp.obj = list(threshold=u[i], estimate=fit2a$par[c(i,i+d)]))
}
#Marginal estimates for scaled negative Dirichlet
if(plots){
  tikz("figure/qqplots.tex",standAlone=TRUE, width=7, height=3.25)
  par(mfrow=c(1,3))
  for(i in 1:d){
    gp.qqplot(x=x[,i], gp.obj = list(threshold=u[i], estimate=fit2b$par[c(i,i+d)]),
              main=paste0("QQ-plot, ",colnames(dat[,i])))
  }
  dev.off()
}

#Alternative fit using the gradient score

#Transform observations to the unit Generalized Pareto scale
pareto.dat <- 1/(1-apply(dat,2, rank, ties.method = "random")/(nrow(dat)+1))
pareto.dat <- matrix(unlist(pareto.dat), ncol=3)
#This does not deal with the matter of clustering, but is used as an illustrative purpose
dir.score.par <- fscore(start=unlist(fit4a$estimate)[7:10], dat=pareto.dat, model="dir", p=20, qu=0.9)$par
negdir.score.par <- fscore(start=unlist(fit4b$estimate)[7:10], dat=pareto.dat, model="negdir", p=20, qu=0.9)$par
#We can compare with the estimates
rbind(negdir.score.par, c(fit4b$estimate$alpha, fit4b$estimate$rho))
#Estimates are comparable, uncertainty diagnostics ignored
#Note: same trick could apply, but unclear what happens given that a bootstrap scheme would
#consist in taking all the data, with more variability.

###############################################################
#### Figures from Extremal Attractors of Liouville Copulas ####
###############################################################

#Some 2d and 3d plots of the angular measure for selected parameter values
plotdens <- function(simplexfn, main, res=520, isLog=TRUE, ...){
  #From Xi'an Og
  # density on a grid with NAs outside, as in image()
  # oldpar <- par(no.readonly = TRUE)
  gride=matrix(NA,ncol=res,nrow=res)
  ww3=ww2=seq(0.01,0.99,le=res)
  for (i in 1:res){
    cur=ww2[i];op=1-cur
    for (j in 1:(res+1-i)){
      gride[i,j]=simplexfn(cbind(cur,ww3[j],op-ww3[j]),...)
      if(!isLog){
        gride[i,j] <- exp(gride[i,j])
      }
      #gride[i,j]=simplexfn(cbind(cur,ww3[j],op-ww3[j]),param=c(1.18,1.19,1.2,1.1),d=3,transform=FALSE)
    }
  }
  subset=(1:length(gride))[is.finite(gride)]
  logride=gride[subset]
  grida=(logride-min(logride))/diff(range(logride))
  grolor=grey.colors(250,start=0,end=0.7)[1+trunc(as.vector(grida)*250)]
  iis=(subset%%res)+res*(subset==res)
  jis=(subset%/%res)+1
  x0 <- (ww3[jis]-ww2[iis])/sqrt(2)
  y0 <- (1-ww3[jis]-ww2[iis])/sqrt(2/3)
  
  #par(mfrow=c(1,2))
  gride[!is.finite(gride)] <- NA
  #image(ww3,ww2,gride,bty="l",xlab="$w_1$",ylab="$w_2$", col=rainbow(250,start=0,end=0.7))
  # preparing the graph
  # setting the limits
  #plot(c(0,1),col="white",axes=FALSE, xlab="",ylab="",#pty="s",
  #     xlim=c(-1.25,1.25)*1.1/sqrt(2),ylim=c(-.2,sqrt(3/2))*1.2, main=main)
  fields::quilt.plot(x0,y0,as.vector(grida),nlevel=250,axes=FALSE,add.legend=FALSE,xlim=c(min(x0)-0.4,max(x0)+0.4),ylim=c(min(y0)-0.2,max(y0)+0.2), main=main)
  polygon(x=c(min(x0)-0.1,max(x0)+0.1,max(x0)+0.1,min(x0)-0.1),y=rep(c(0.005,0-0.2),each=2),col="white",border =NA)
  polygon(x=c(min(x0),rep(0,2),rep(min(x0)-0.2,3)),y=c(min(y0),max(y0),rep(max(y0)+0.2,2),max(y0),min(y0)),col='white',border =NA)
  polygon(x=c(max(x0),rep(0,2),rep(max(x0)+0.2,3)),y=c(min(y0),max(y0),rep(max(y0)+0.2,2),max(y0),min(y0)),col='white',border = NA)
  #segments(x0=min(x0),x1=max(x0),y0=min(y0),col=1,lwd=2)
  #segments(x0=min(x0),x1=0,y0=min(y0),y1=max(y0),col=1,lwd=2)
  #segments(x0=0,x1=max(x0),y0=max(y0),y1=min(y0),col=1,lwd=2)
  text(x=min(x0)-0.15, y=-0.1, "$w_1=1$",pos = 4)
  text(x=max(x0)+0.15, y=-0.1, "$w_2=1$",pos = 2)
  text(x=0, y=max(y0)+0.15, "$w_3=1$", pos=1)
  #par(oldpar)
}

if(plots){
  tikz("figure/angdens.tex",standAlone=TRUE, width=7, height=2.8)
  par(mfrow=c(1,3), mar=c(1,0.1,3,0.1))
  plotdens(specdens,param=c(fit5a$est$alpha, fit5a$est$rho),d=3,transform=FALSE, isLog=FALSE, main="Negative logistic",res=300, model="dir")
  plotdens(specdens,param=c(fit5b$est$alpha, fit5b$est$rho),d=3,transform=FALSE, isLog=FALSE, main="Logistic",res=300, model="negdir")
  plotdens(specdens,param=c(fit2b$est$alpha, fit2b$est$rho),d=3,transform=FALSE, isLog=FALSE, main="scaled\n Dirichlet",res=300, model="negdir")
  dev.off()
}


#Tables 1 and 2
library(xtable)
se <- c(sapply(1:d, function(i){mev::gp.fit(x[,i],u[i])$std.err}))
se <- se[c(1,3,5,2,4,6)]
fit5b$se


table.res <- matrix(nrow=5, ncol=6)
table.res[1,1:3] <- sapply(1:3, function(ind){paste0(round(fit2b$par[ind],1), " (",round(fit2b$se[ind],1),")")})
table.res[2,1:3] <- sapply(1:3, function(ind){paste0(round(fit5a$par[ind],1), " (",round(fit5a$se[ind],1),")")})
table.res[3,1:3] <- sapply(1:3, function(ind){paste0(round(fit5b$par[ind],1), " (",round(fit5b$se[ind],1),")")})
table.res[4,1:3] <- sapply(1:3, function(ind){paste0(round(fit1a$par[ind],1), " (",round(fit1a$se[ind],1),")")})
table.res[5,1:3] <- sapply(1:3, function(ind){paste0(round(c(scale, shape)[ind],1), " (",round(se[ind],1),")")})

table.res[1,4:6] <- sapply(4:6, function(ind){paste0(round(fit2b$par[ind],2), " (",round(fit2b$se[ind],2),")")})
table.res[2,4:6] <- sapply(4:6, function(ind){paste0(round(fit5a$par[ind],2), " (",round(fit5a$se[ind],2),")")})
table.res[3,4:6] <- sapply(4:6, function(ind){paste0(round(fit5b$par[ind],2), " (",round(fit5b$se[ind],2),")")})
table.res[4,4:6] <- sapply(4:6, function(ind){paste0(round(fit1a$par[ind],2), " (",round(fit1a$se[ind],2),")")})
table.res[5,4:6] <- sapply(4:6, function(ind){paste0(round(c(scale, shape)[ind],2), " (",round(se[ind],2),")")})


table.res2 <- matrix(nrow=5, ncol=4)
table.res2[1,-4] <- sapply(7:9, function(ind){paste0(round(fit2b$par[ind],2), " (",round(fit2b$se[ind],2),")")})
table.res2[1,4] <- sapply(10, function(ind){paste0("$-$",round(fit2b$par[ind],2), " (",round(fit2b$se[ind],2),")")})
table.res2[2,4] <- sapply(7, function(ind){paste0(round(fit5a$par[ind],2), " (",round(fit5a$se[ind],2),")")})
table.res2[3,4] <- sapply(7, function(ind){paste0(round(fit5b$par[ind],2), " (",round(fit5b$se[ind],2),")")})
table.res2[4,-4] <- sapply(7:9, function(ind){paste0(round(fit1a$par[ind],2), " (",round(fit1a$se[ind],2),")")})
table.res2[5,1:3] <- sapply(1:3, function(ind){paste0(round(negdir.score.par[ind],2))})
table.res2[5,4] <- sapply(4, function(ind){paste0("$-$",round(negdir.score.par[ind],2))})
table.res2[2:3,-4] <- "1"
table.res2[4,4] <- "1"

rownames(table.res) <- c("Scaled Dirichlet","Neg. logistic","Logistic","ext. Dirichlet","Marginal")
colnames(table.res) <- c(sapply(1:3, function(i){paste0("$\\eta_",i,"$")}), sapply(1:3, function(i){paste0("$\\xi_",i,"$")}))

rownames(table.res2) <- c("Scaled Dirichlet","Neg. logistic","Logistic","ext. Dirichlet","Gradient score")
colnames(table.res2) <- c(sapply(1:3, function(i){paste0("$\\alpha_",i,"$")}), "$\\rho$")
library(xtable)
#Removing the scale parameters
tab1 <- xtable(table.res,caption = "Generalized Pareto parameter estimates and standard errors (in parenthesis) for the trivariate river example for four different models. ",label="tab:margpar")
tab2 <- xtable(table.res2,caption = "Dependence parameters estimates and standard errors (in parenthesis) for the trivariate river example.",label="tab:deppar")
#The three first estimates were obtained using the routine in ``fmvcpot'' and the last via ``fscore'', with standard errors estimated via a nonparametric boostrap.
print(tab1, file="Table1.tex", type="latex", sanitize.text.function=function(str){str}, math.style.negative=TRUE, booktabs=TRUE)
print(tab2, file="Table2.tex", type="latex", sanitize.text.function=function(str){str}, math.style.negative=TRUE, booktabs=TRUE)


###Figures 1, 2 and 3

#Some 2d and 3d plots of the angular measure for selected parameter values
plotdens <- function(simplexfn, main, res=520, isLog=TRUE, ...){
  #From Xi'an Og
  # density on a grid with NAs outside, as in image()
  # oldpar <- par(no.readonly = TRUE)
  gride=matrix(NA,ncol=res,nrow=res)
  ww3=ww2=seq(0.01,0.99,le=res)
  for (i in 1:res){
    cur=ww2[i];op=1-cur
    for (j in 1:(res+1-i)){
      gride[i,j]=simplexfn(cbind(cur,ww3[j],op-ww3[j]),...)
      if(!isLog){
        gride[i,j] <- exp(gride[i,j])
      }
    }
  }
  subset=(1:length(gride))[is.finite(gride)]
  logride=gride[subset]
  grida=(logride-min(logride))/diff(range(logride))
  grolor=grey.colors(250,start=0,end=0.7)[1+trunc(as.vector(grida)*250)]
  iis=(subset%%res)+res*(subset==res)
  jis=(subset%/%res)+1
  x0 <- (ww3[jis]-ww2[iis])/sqrt(2)
  y0 <- (1-ww3[jis]-ww2[iis])/sqrt(2/3)


  gride[!is.finite(gride)] <- NA
  fields::quilt.plot(x0,y0,as.vector(grida),nlevel=250,axes=FALSE,add.legend=FALSE,xlim=c(min(x0)-0.4,max(x0)+0.4),ylim=c(min(y0)-0.2,max(y0)+0.2), main=main)
  polygon(x=c(min(x0)-0.1,max(x0)+0.1,max(x0)+0.1,min(x0)-0.1),y=rep(c(0.005,0-0.2),each=2),col="white",border =NA)
  polygon(x=c(min(x0),rep(0,2),rep(min(x0)-0.2,3)),y=c(min(y0),max(y0),rep(max(y0)+0.2,2),max(y0),min(y0)),col='white',border =NA)
  polygon(x=c(max(x0),rep(0,2),rep(max(x0)+0.2,3)),y=c(min(y0),max(y0),rep(max(y0)+0.2,2),max(y0),min(y0)),col='white',border = NA)
  text(x=min(x0)-0.15, y=-0.1, "$w_1=1$",pos = 4)
  text(x=max(x0)+0.15, y=-0.1, "$w_2=1$",pos = 2)
  text(x=0, y=max(y0)+0.15, "$w_3=1$", pos=1)
}
if(plots){

  #Trivariate angular density
  tikz("spec_dens_3d.tex",width=7, height=6.5, standAlone=TRUE)
  par(mfrow=c(2,2), mar=c(1,0.1,3,0.1))
  plotdens(specdens,param=c(1,0.5,0.2,0.2),d=3,transform=FALSE, isLog=FALSE, main="scaled\n Dirichlet",res=200, model="dir")
  plotdens(specdens,param=c(1.25,2,1,0.4),d=3,transform=FALSE, isLog=FALSE, main="scaled negative\nDirichlet",res=200, model="negdir")
  plotdens(specdens,param=c(0.2,0.2,0.2,0.2),d=3,transform=FALSE, isLog=FALSE, main="scaled\n Dirichlet",res=200, model="dir")
  plotdens(specdens,param=c(1.25,1.25,1.25,0.4),d=3,transform=FALSE, isLog=FALSE, main="scaled negative\nDirichlet",res=200, model="negdir")
  dev.off()
  graphics.off()

  #Bivariate angular density

  tikz("spec_dens_2d.tex",width=7, height=3.5, standAlone=TRUE)
  #Bivariate spectral density plots
  xseq <- seq(0,1,by=0.001)
  par(mfrow=c(1,2))
  plot(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="dir", param=c(2,0.5,0.8),d=2,transform=FALSE)})),
    type="l",xlab="$w$",ylab="Angular density",bty="l")
  lines(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="dir", param=c(0.1,0.1,0.25),d=2,transform=FALSE)})),
    col=2, lty=2,lwd=1.5)
  lines(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="dir", param=c(0.5,0.5,0.25),d=2,transform=FALSE)})),
    col=4, lty=3,lwd=2)
  title("scaled Dirichlet")
  #Could have used `curve' here too
  plot(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="negdir",param=c(2,0.5,0.25),d=2,transform=FALSE)})),
    type="l",xlab="$w$",ylab="Angular density",bty="l")
  lines(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="negdir",param=c(0.4,0.4,0.25),d=2,transform=FALSE)})),
    col=2, lty=2,lwd=1.5)
  lines(xseq,exp(sapply(xseq, function(x){
    specdens(dat=cbind(x,1-x), model="negdir", param=c(0.5,0.5,0.25),d=2,transform=FALSE)})),
    col=4, lty=3,lwd=2)

  title("scaled negative Dirichlet")

  dev.off()


  #Pickands dependence function

  tikz("pickands_plots.tex",standAlone=TRUE,width=7,height=3.5)
  par(mfrow=c(1,2),mai = c(1, 0.7, 0.7, 0.1))
  lcopula::pickands.plot(alpha = c(2,0.5),rho=0.8, plot.new=TRUE,CDA="S",tikz=TRUE)
  lcopula::pickands.plot(alpha = c(0.1,0.1),rho=0.25, plot.new=FALSE, CDA="S", tikz=TRUE, col=2, lwd=1.5, lty=2)
  lcopula::pickands.plot(alpha = c(0.5,0.5),rho=0.25, plot.new=FALSE, CDA="S", tikz=TRUE, col=4, lwd=2, lty=3)
  title("scaled Dirichlet")
  lcopula::pickands.plot(alpha = c(2,0.5),rho=0.25, plot.new=TRUE,CDA="C",tikz=TRUE)
  lcopula::pickands.plot(alpha = c(0.4,0.4),rho=0.25, plot.new=FALSE, CDA="C", tikz=TRUE, col=2, lwd=1.5, lty=2)
  lcopula::pickands.plot(alpha = c(0.5,0.5),rho=0.25, plot.new=FALSE, CDA="C", tikz=TRUE, col=4, lwd=2, lty=3)
  title("scaled negative Dirichlet")
  dev.off()
}
