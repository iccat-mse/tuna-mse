library(FLBRP)
library(FLash)
library(biodyn)

library(RSQLite)
library(DBI)

dirDat='~/Desktop/MEGA/papers/submitted/tuna-mse/data'
dirRes='~/Desktop/MEGA/papers/submitted/tuna-mse/resubmission/data'

load(paste(dirRes,"/mou.RData",     sep=""))
load(paste(dirRes,"/srDev.RData",   sep=""))
load(paste(dirRes,"/options.RData", sep=""))
load(paste(dirDat,"/BRPs.RData",    sep=""))
load(paste(dirDat,"/priors.RData",  sep=""))
load(paste(dirRes,"/newScen.RData", sep=""))

start=50; end=100; interval =3
nits =100 

trendQ=FLQuant(1,dimnames=list(year=1:end))
omega =1
uCV   =0.3
db    =NULL

# dataPoorMSE=function(om,br,srDev,priors,
#                      start=50, end =100,interval=3,
#                      nits =100,seed=7890,
#                      omega =1,
#                      trendQ=FLQuant(1,dimnames=list(year=1:end)),
#                      uCV=0.3,
#                      db=NULL){
  
  set.seed(seed)
      
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

