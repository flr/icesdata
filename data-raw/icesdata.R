load("C:/active/flr/icesdata/data/raw/Updated_stks_n82_R0_updated202602_15.Rdata")

icesdata=ICESStocks

rm(ICESStocks)

for (.id in subset(info,Model=="SMS")$.id)
  icesdata[[.id]]=qapply(icesdata[[.id]], function(x) {
    dimnames(x)[3:5]=dimnames(stock.n(icesdata[[.id]]))[3:5]
    
    dimnames(x)[[3]]="unique"
    dimnames(x)[[4]]="all"
    dimnames(x)[[5]]="unique"
    x})

for(.id in names(icesdata)){
  small=min(catch.n(icesdata[[.id]]),na.rm=TRUE)*1e-6
  catch.n(   icesdata[[.id]])[is.na(catch.n(   icesdata[[.id]]))]=small
  landings.n(icesdata[[.id]])[is.na(landings.n(icesdata[[.id]]))]=small
  discards.n(icesdata[[.id]])[is.na(catch.n(   icesdata[[.id]]))]=0.0}

for (san in c("san.sa.1r","san.sa.2r","san.sa.3r","san.sa.4"))
  icesdata[[san]]=qapply(icesdata[[san]], function(x) {
    dimnames(x)$unit="unique"
    x})

save(icesdata,file="C:/active/flr/icesdata/data/icesdata.Rdata")

ctc1903=read.csv("C:/active/flr/icesdata/data-raw/inputs/NorthSea_stocks.csv")[,1:12] 
ctc1903$Discards[is.na(ctc1903$Discards)]=0
names(ctc1903)[c(3:4,6)]=c(".id","year","catch")
names(ctc1903)=tolower(names(ctc1903))
save(ctc1903,file="C:/active/flr/icesdata/data/ctc1903")
