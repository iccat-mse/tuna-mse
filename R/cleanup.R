library(FLCore)
library(FLBRP)
library(FLash)
library(biodyn)
library(reshape)

library(RSQLite)
library(DBI)

drv=dbDriver("SQLite")

fls=c("non21.db","rerun.db","bndf.db","bndTac.db")
dir   ="/home/laurie/MEGAsync/papers/-rfmo-mse"
dirOld="/home/laurie/Desktop/MEGA/papers/submitted/tuna-mse/resubmission"

con=mlply(fls, function(db)  dbConnect(drv, dbname=file.path(dirOld,"data",db)))
om =llply(con, dbReadTable, "om")
mou=llply(con, dbReadTable, "mou")
mp =llply(con, dbReadTable, "mp")

om[[ 2]]=subset(om[[2]], OM==21)
mou[[2]]=subset(mou[[2]],OM==21)
mp[[ 2]]=subset(mp[[2]], OM==21)
om[[ 2]]=om[[ 2]][!duplicated(om[[ 2]][,1:7]),]
mou[[2]]=mou[[2]][!duplicated(mou[[2]][,1:7]),]
mp[[ 2]]=mp[[ 2]][!duplicated(mp[[ 2]][,c("OM","SA","HCR","omega","qTrend","p","year","iter")]),]

om =rbind(cbind(rbind(om[[1]],om[[2]]),   boundF=FALSE,boundTac=FALSE),
          cbind(      om[[3]],            boundF=TRUE, boundTac=FALSE),
          cbind(      om[[4]],            boundF=FALSE,boundTac=TRUE))
mp =rbind(cbind(rbind(mp[[1]],mp[[2]]),   boundF=FALSE,boundTac=FALSE),
          cbind(mp[[3]],                  boundF=TRUE, boundTac=FALSE),
          cbind(mp[[4]],                  boundF=FALSE,boundTac=TRUE))
mou=rbind(cbind(rbind(mou[[1]],mou[[2]]), boundF=FALSE,boundTac=FALSE),
          cbind(mou[[3]],                 boundF=TRUE, boundTac=FALSE),
          cbind(mou[[4]],                 boundF=FALSE,boundTac=TRUE))
save(om, file=file.path(dir,"data","om.RData"), compress="xz")
save(mou,file=file.path(dir,"data","mou.RData"),compress="xz")
save(mp, file=file.path(dir,"data","mp.RData"), compress="xz")

tmp=melt(with(om,table(OM,SA,HCR,bound,qTrend,omega)))

ggplot(subset(om,year>55))+
  geom_boxplot(aes(x=factor(year),y=harvest))+
  facet_grid(HCR~SA)

ggplot(subset(om2,year>55&OM==5))+
  geom_boxplot(aes(x=factor(year),y=harvest))+
  facet_grid(HCR~SA)

ggplot(bc)+
  geom_boxplot(aes(x=factor(year),y=ssb))
