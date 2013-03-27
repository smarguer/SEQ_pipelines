toLINEAR <- function(f1){
out=2^f1
out[which(out==0.1)]=0
return (out)
}

toLOG <- function(f1){
f1[which(f1==0)]=0.1
out=log2(f1)
return(out)
}

toLINEAR_norm <- function(f1){
out[which(out==0)]=1
out=2^f1
out[which(out==2)]=0
return (out)
}

toLOG_norm<- function(f1){
f1[which(f1==0)]=2
out=log2(f1)
f1[which(f1==1)]=0

return(out)
}


