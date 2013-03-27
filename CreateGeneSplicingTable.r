CreateGeneSplicingTable=function(tr,ei)
{
 out=data.frame()
 for (i in 1:nrow(ei))
 {
  out[i,1]=ei[i,1]
  out[i,2]=ei[i,6]
  out[i,3]=sum(tr[grep(ei[i,1], tr[,1]),7], na.rm=TRUE)
  out[i,4]=sum(tr[grep(ei[i,1], tr[,1]),8], na.rm=TRUE)
  out[i,5]=ei[i,7]
  out[i,6]=ei[i,8]
  out[i,7]=ei[i,11]
  out[i,8]=ei[i,12]
  out[i,9]=round((out[i,3]/out[i,4]),2)
  out[i,10]=round((out[i,3]/(out[i,3]+out[i,4])*100),2)
  out[i,11]=round((out[i,5]/out[i,6]),2)
 }
 colnames(out)=c(colnames(ei)[c(1,6)],"tr","int","signal.exonic","signal.intronic","RC.exonic","RC.intronic","sp.eff","perc.sp","ei.ratio")
 out=out[which(is.na(out[,"signal.intronic"])==FALSE),]
 return(out)
}



