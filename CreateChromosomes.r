CreateChromosomes<-function(f1.dir,f1,depth)
{

#source("format.bases.r")
#source("toLINEAR.r")

#cat("******")
#cat("\n")

chr1_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr1.bases",sep=""))
#cat("*")
chr2_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr2.bases",sep=""))
#cat("*")
chr3_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr3.bases",sep=""))
#cat("*")
chr1_w5t=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLtrs.GEN.col.chr1.bases",sep=""))
#cat("*")
chr2_w5t=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLtrs.GEN.col.chr2.bases",sep=""))
#cat("*")
chr3_w5t=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLtrs.GEN.col.chr3.bases",sep=""))
#cat("*")
#cat("\n")

chr1=cbind(toLOG(chr1_w5),toLOG(chr1_w5t))
chr2=cbind(toLOG(chr2_w5),toLOG(chr2_w5t))
chr3=cbind(toLOG(chr3_w5),toLOG(chr3_w5t))
save(list=c("chr1","chr2","chr3"), file=paste(f1,".EACH.rda",sep=""))
rm(chr1,chr2,chr3)

chr1=chr1_w5+chr1_w5t
chr2=chr2_w5+chr2_w5t
chr3=chr3_w5+chr3_w5t
chr1=toLOG(chr1)
chr2=toLOG(chr2)
chr3=toLOG(chr3)
save(list=c("chr1","chr2","chr3"), file=paste(f1,".rda",sep=""))

chr1=toLINEAR(chr1)
chr2=toLINEAR(chr2)
chr3=toLINEAR(chr3)

write.table(chr1, file=paste(f1,".chr1",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(chr2, file=paste(f1,".chr2",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(chr3, file=paste(f1,".chr3",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)

chr1=chr1/depth
chr2=chr2/depth
chr3=chr3/depth

chr1=toLOG(f1=chr1)
chr2=toLOG(f1=chr2)
chr3=toLOG(f1=chr3)
save(list=c("chr1","chr2","chr3"), file=paste("NORM_",f1,".rda",sep=""))

rm(list=ls(pattern="chr"))
rm(depth)

}
CreateChromosomes.CHIP<-function(f1.dir,f1,depth)
{

#source("format.bases.r")
#source("toLINEAR.r")

#cat("******")
#cat("\n")

chr1_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr1.bases",sep=""))
#cat("*")
chr2_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr2.bases",sep=""))
#cat("*")
chr3_w5=format.bases.nolog(f0=paste(f1.dir,f1,".5mis.ALLchr.col.chr3.bases",sep=""))
#cat("*")
chr1_w5[1:(length(chr1_w5)/2)]=chr1_w5[1:(length(chr1_w5)/2)]+chr1_w5[((length(chr1_w5)/2)+1):length(chr1_w5)]
chr2_w5[1:(length(chr2_w5)/2)]=chr2_w5[1:(length(chr2_w5)/2)]+chr2_w5[((length(chr2_w5)/2)+1):length(chr2_w5)]
chr3_w5[1:(length(chr3_w5)/2)]=chr3_w5[1:(length(chr3_w5)/2)]+chr3_w5[((length(chr3_w5)/2)+1):length(chr3_w5)]
chr1_w5[((length(chr1_w5)/2)+1):length(chr1_w5)]=chr1_w5[1:(length(chr1_w5)/2)]
chr2_w5[((length(chr2_w5)/2)+1):length(chr2_w5)]=chr2_w5[1:(length(chr2_w5)/2)]
chr3_w5[((length(chr3_w5)/2)+1):length(chr3_w5)]=chr3_w5[1:(length(chr3_w5)/2)]


chr1=toLOG(chr1_w5)
chr2=toLOG(chr2_w5)
chr3=toLOG(chr3_w5)
save(list=c("chr1","chr2","chr3"), file=paste(f1,".rda",sep=""))

chr1=toLINEAR(chr1)
chr2=toLINEAR(chr2)
chr3=toLINEAR(chr3)

write.table(chr1, file=paste(f1,".chr1",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(chr2, file=paste(f1,".chr2",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(chr3, file=paste(f1,".chr3",sep=""), row.names=FALSE,col.names=FALSE,quote=FALSE)

chr1=chr1/depth
chr2=chr2/depth
chr3=chr3/depth

chr1=toLOG(f1=chr1)
chr2=toLOG(f1=chr2)
chr3=toLOG(f1=chr3)
save(list=c("chr1","chr2","chr3"), file=paste("NORM_",f1,".rda",sep=""))

rm(list=ls(pattern="chr"))
rm(depth)

}
