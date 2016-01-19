library(FLCore)
library(FLBRP)
library(FLash)
library(biodyn)
library(ggplotFL)

library(RSQLite)
library(DBI)

dirDat='~/Desktop/MEGA/papers/submitted/tuna-mse/data'
db    =file.path(dirDat,"mseTest.db")

##Data sets
# The design grid for all interactions
load(paste(dirDat,"/design.RData",   sep=""))

# FLBRPs and FLStocks based on design
load(paste(dirDat,"/OMs.RData",      sep=""))
load(paste(dirDat,"/BRPs.RData",     sep=""))

# Priors for biodyn based on design
load(paste(dirDat,"/priors.RData",   sep=""))
priors=priors[, dimnames(biodyn()@priors)$params]

## OM Scenarios
options=list(OM =expand.grid(OM=seq(128),ar=c(0.0,0.5)),
             HCR=expand.grid(ftar=c(.5,.75),blim=c(.3,.4),btrig=c(0.6,0.8)),
             SA =expand.grid(p=c(1,2),prr=c("None","r","k","fmsy","bmsy"),stringsAsFactors=FALSE),
             OEM=expand.grid(omega=c(1,.75),qTrend=c(0,.02)))


## Range and iters
start=50; end=1500; interval=3; rcvPeriod=3
nits =100; seed=7890 

set.seed(seed)
srDev=FLQuants("0"=exp(lh:::noise(nits, FLQuant(0,dimnames=list(year=1:(end+interval))),
                         b=0,sd=.3)))

set.seed(seed)
srDev[["0.5"]]   =exp(lh:::noise(nits, FLQuant(0,dimnames=list(year=1:(end+interval))),
                         b=0.3,sd=0.5))

uDev   =rlnorm(nits,FLQuant(0,dimnames=list(year=1:end)),0.3)

scen=do.call("expand.grid",llply(options,function(x) seq(length(dimnames(x)[[1]]))))
names(scen)=paste("i",names(scen),sep="")

m_ply(subset(scen,(iOM%in%c(21,5,17,22,23,29,53,85)[1]))[1,], 
      function(iOM,iSA,iHCR,iOEM){

  eql=BRPs[[options$OM[iOM,"OM"]]]
  srD=srDev[[ac(options$OM[iOM,"ar"])]]

  #### OM
  om=OMs[[options$OM[iOM,"OM"]]]
  om=fwdWindow(om,end=end+interval,eql)
    
  ## F in recovery period
  rcv.=seq(c(fbar(om)[,ac(start)]),c(FLBRP:::refpts(eql)["msy","harvest"])*options$HCR[iHCR,"ftar"],length.out=rcvPeriod+1)
  rcv =FLQuant(rcv.,dimnames=dimnames(fbar(om)[,ac(start+0:rcvPeriod)]))
  om  =FLash:::fwd(om,f=rcv, sr=eql)
    
  ## F in longterm
  lgt =propagate(FLQuant(rcv.[rcvPeriod+1],dimnames=list(year=(start+rcvPeriod+1):(end+interval))),nits)
  om  =propagate(om,nits)
  om  =FLash:::fwd(om,f=lgt, sr=eql)
    
  ## Add stochastcity
  om =FLash:::fwd(om,f =rlnorm(nits,log(fbar(FLCore:::iter(om,1))[,ac(2:(end+interval))]),0.1),
                     sr=eql,sr.residuals=srD)

  ## save projection for comparison later
  mou=om
  trendQ=FLQuant(c(rep(1,start-1),cumprod(rep(1+options$OEM[iOEM,"qTrend"],end-start+1))),
                 dimnames=list(year=1:end))

  #### MP
  ## SA Control
  ctrl=with(priors[iOM,],biodyn:::controlFn(r=r,k=k,p=p,b0=b0))
  if (options$SA[ iSA, "p"]==1) ctrl["p","val"]=1
  
  ## SA Priors
  prArg=data.frame(a=unlist(c(priors[iOM,])),b=unlist(c(priors[iOM,]))*0.5,weight=1)
  prArg=alply(prArg,1,cbind)
  names(prArg)=names(priors[iOM,])
  prs=do.call(priorFn,prArg[options$SA[iSA,2]])
  
  #### MSE
  res=do.call(mseBiodyn,
             list(om,eql=eql,srDev=srD,
                  control=ctrl,priors=prs,
                  start=start+rcvPeriod,end=end,interval=3,
                  ftar=options$HCR[iHCR,"ftar"],
                  blim=options$HCR[iHCR,"blim"],
                  btrig=options$HCR[iHCR,"btrig"],
                  uDev  =uDev,
                  qTrend=trendQ,
                  omega =options$OEM[iOEM,"omega"],
                  maxF  =1,
                  refB=FLBRP:::refpts(eql)["msy","biomass"]))
  
  #### Save results
  drv=dbDriver("SQLite")
  con=dbConnect(drv, dbname=db)
    
  mou=cbind(OM=iOM,SA=iSA,HCR=iHCR,OEM=iOEM,biodyn:::tseries(res$mou,eql))
  om =cbind(OM=iOM,SA=iSA,HCR=iHCR,OEM=iOEM,biodyn:::tseries(res$om, eql))
  mp =cbind(OM=iOM,SA=iSA,HCR=iHCR,OEM=iOEM,res$mp)
    
  dbWriteTable(con, "mou", mou,append=TRUE)
  dbWriteTable(con, "om",  om, append=TRUE)
  dbWriteTable(con, "mp",  mp, append=TRUE)
  
  dbDisconnect(con)})

drv=dbDriver("SQLite")
con=dbConnect(drv, dbname=db)
om =dbReadTable(con, "om")
