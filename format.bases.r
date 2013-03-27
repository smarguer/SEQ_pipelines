format.bases<-function(f0)
{
test=read.delim(f0, header=F)
test=test[order(test[,1]),]
test=test[,2]
test[which(test==0)]=0.1
test=log2(test)
return(test)
}

format.bases.nolog<-function(f0)
{
test=read.delim(f0, header=F)
test=test[order(test[,1]),]
test=test[,2]
#test[which(test==0)]=0.1
#test=log2(test)
return(test)
}

