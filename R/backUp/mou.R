##sets up MOUs, i.e. projections without feedback

library(FLBRP)
library(FLash)

dirDat='~/Desktop/MEGA/papers/submitted/tuna-mse/data'
dirRes='~/Desktop/MEGA/papers/submitted/tuna-mse/resubmission/data'

load(paste(dirDat,"/OMs.RData",    sep=""))
load(paste(dirDat,"/BRPs.RData",   sep=""))
load(paste(dirDat,"/design.RData", sep=""))
load(paste(dirDat,"/priors.RData", sep=""))
load(paste(dirDat,"/omKey.RData",  sep=""))

setMou<-function(om,br,srDev,
              ftar =0.75,
              fCV  =0.3,
              start=50,end=103,rcvPeriod=5){
  
  ## Recovery
  rcv =seq(c(fbar(om)[,ac(start)]),c(FLBRP:::refpts(br)["msy","harvest"])*ftar,length.out=rcvPeriod+1)
  rcv =FLQuant(rcv,dimnames=dimnames(fbar(om)[,ac(start+0:rcvPeriod)]))
  om  =FLash:::fwd(om,f=rcv,sr=br)
  
  ## F in longterm
  om  =fwdWindow(om,br,end=end)
  lgt =FLQuant(c(rcv[,rcvPeriod+1]),dimnames=list(year=(start+rcvPeriod+1):(end)))
  om  =FLash:::fwd(om,f=lgt, sr=br)
  
  harvest(om)=rlnorm(dim(srDev)[6],log(harvest(om)),fCV)
  units(harvest(om))="f"
  apex=c(ages(catch.sel(br))[catch.sel(br)==c(fapex(catch.sel(br)))])
  range(om)[c("minfbar","maxfbar")]=apex
  res=fwd(om,f=fbar(om)[,-1],sr=br,sr.residuals=srDev)
  
  res}

## OM
options=list(OM =rbind(data.frame(i=as.numeric(names(omKey)),ar=0.0),
                       data.frame(i=21,                      ar=0.5)))

## HCRs
options$HCR =expand.grid(ftar=c(.5,.75),blim=c(.3,.4),btrig=c(0.6,0.8))

## SA
options$SA =rbind(data.frame(p=1,prr=c("None","r","k","fmsy","bmsy"),stringsAsFactors=FALSE),
                  data.frame(p=2,prr=c("None","r","k","fmsy"),       stringsAsFactors=FALSE))

options$OEM =expand.grid(omega=c(0,.75),qTrend=c(0,.02))

save(options,file=paste(dirRes,"/options.RData", sep=""))

start=50; end=100; interval =3
nits =100 

set.seed(7890)
srDev=FLQuants("0"=exp(lh:::noise(nits, FLQuant(0,dimnames=list(year=1:(end+interval))),
                         b=0,sd=.3)))

set.seed(7890)
srDev[["0.5"]]   =exp(lh:::noise(nits, FLQuant(0,dimnames=list(year=1:(end+interval))),
                         b=0.3,sd=0.5))

##OM with perfect management
mou=FLStocks(mlply(options$OM, function(i,ar)
                setMou(OMs[[i]],BRPs[[i]],srDev[[ac(ar)]])))

save(mou,  file=paste(dirRes,"mou.RData",  sep="/"))
save(srDev,file=paste(dirRes,"srDev.RData",sep="/"))

