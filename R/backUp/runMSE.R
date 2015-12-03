library(FLBRP)
library(FLash)
library(biodyn)
library(FLife)

library(RSQLite)
library(DBI)

source('~/Desktop/MEGA/papers/submitted/tuna-mse/submission/R/biodyn-tseries.R')
source('~/Desktop/MEGA/papers/submitted/tuna-mse/submission/R/rand-noise.R')
source('~/Desktop/MEGA/papers/submitted/tuna-mse/submission/R/biodyn-prior.R')
#source('~/Desktop/MEGA/papers/submitted/tuna-mse/submission/R/biodyn-hcr.R')
source('~/Desktop/MEGA/papers/submitted/tuna-mse/submission/R/biodyn-mse.R')

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

save(options,file=file.path(dirRes,"options.RData"))

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

save(mou,file=paste(dirRes,"mou.RData",sep="/"))

trendQ=FLQuant(1,dimnames=list(year=1:end))
omega =1
uCV   =0.3
db    =NULL

#### Set options
prr  =options$SA[ iSA, "prr"]
ftar =options$HCR[iHCR,"ftar"]
btrig=options$HCR[iHCR,"btrig"]
blim =options$HCR[iHCR,"blim"]

# dataPoorMSE=function(om,br,srDev,priors,
#                      start=50, end =100,interval=3,
#                      nits =100,seed=7890,
#                      omega =1,
#                      trendQ=FLQuant(1,dimnames=list(year=1:end)),
#                      uCV=0.3,
#                      db=NULL){
  
  set.seed(seed)
  
  #### OM
  if (dims(om)$maxyear < end+interval)
    om =fwdWindow(om,end=end+interval,br)
  
  #### MP
  ## SA Options
  ctrl=with(priors[iOM,c("r","k","p","b0")],biodyn:::controlFn(r=r,k=k,p=p,b0=b0))
  if (options$SA[ iSA, "p"]==1) ctrl["p","val"]=1
    
  val=priors[options$OM[iOM,"i"],prr]
  prArg=list(c(weight=1,a=val,b=val*0.5))
  names(prArg)=prr
  prrs=priorFn(prArg)
  
  #### MSE
  res=mseBiodyn(om,br,srDev,
                control=ctrl,priors=prrs,
                start=start+rcvPeriod,end=end,interval=3,
                ftar=ftar,blim=blim,btrig=btrig,
                uDev   =0.3,trendQ=trendQ,
                omega =omega,
                refB  =FLBRP:::refpts(br)["msy","biomass"])

  #### Save results
  if (is.null(db)) return(res)
  
  drv=dbDriver("SQLite")
  con=dbConnect(drv, dbname=db)
  
  mou=cbind(OM=iOM,SA=iSA,HCR=iHCR,tseries(res$mou,br))
  om =cbind(OM=iOM,SA=iSA,HCR=iHCR,tseries(res$om,br))
  mp =cbind(OM=iOM,SA=iSA,HCR=iHCR,res$mp)
  #oem=cbind(OM=iOM,SA=iSA,HCR=iHCR,model.frame(res$oem,drop=T))
  
  dbWriteTable(con, "mou", mou,append=TRUE)
  dbWriteTable(con, "om",  om, append=TRUE)
  dbWriteTable(con, "mp",  mp, append=TRUE)
  #dbWriteTable(con, "oem", oem,append=TRUE)
  
  dbDisconnect(con)}

trendQ=1+FLQuant(cumsum(rep(0.01,end)),
                  dimnames=list(year=1:end))

for (i in 21){
  res=m_ply(expand.grid(iOM=i,iSA =c(1:2,6:7)[4],
                              iHCR=c(1,2,7,8)[3]), 
      function(iOM,iSA,iHCR,options,om,br,priors,db) {
           cat(paste("  OM\t",iOM,":\tSA\t",iSA,":\tHCR\t",iHCR,
           "\n  ======================================",sep=""))
        dataPoorMSE(i,iSA,iHCR,options,om,br,priors,
                    db=NULL,
                    trendQ=trendQ,
                    omega=1)},
            options=options,
            om     =OMs[[      options$OM[i,"i" ]]],
            br     =BRPs[[     options$OM[i,"i" ]]],
            priors=priors,
            db=paste(dirDat,"trendQ",sep="/"))}

db="/home/laurie/Desktop/MEGA/papers/submitted/tuna-mse/data/breakIt6"
db =paste(dirDat,"testOem",sep="/")
drv=dbDriver("SQLite")
con=dbConnect(drv, dbname=db)
oem =dbGetQuery(con,"select * from oem")


oem2(om,cv,trendQ,omega,refB=FLBRP:::refpts(br)["msy","biomass"],
         fishDepend=TRUE)

