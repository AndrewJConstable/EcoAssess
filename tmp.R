outName<-"/perm_storage/home/acon/MEASO/results.rds"
tmp<-c(1:10)
resPeriod<-tmp*100
saveRDS(resPeriod, outName)

#res<-readRDS(outName)