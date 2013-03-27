computeRPKM<-function(f0, gff=gff, gff.name, depth)
{
 j=0
 Chr=read.delim(paste(f0,"ALLchr","col","mapIV",gff.name,"READSFEATURETABLE",sep='.'), stringsAsFactors=F)
 Trs=read.delim(paste(f0,"ALLtrs","GEN","col","mapV",gff.name,"READSFEATURETABLE",sep='.'), stringsAsFactors=F)
 gChr=read.delim(paste(f0,"ALLchr","col","mapIV",gff.name,"READSGENETABLE",sep='.'), stringsAsFactors=F)
 gTrs=read.delim(paste(f0,"ALLtrs","GEN","col","mapV",gff.name,"READSGENETABLE",sep='.'), stringsAsFactors=F)
 colnames(gChr)=c("Name","chr","strand","ORF.start","ORF.end","reads","start","end","AS.reads","AS.start","AS.end")
 colnames(gTrs)=c("Name","chr","strand","ORF.start","ORF.end","reads","start","end","AS.reads","AS.start","AS.end")
 gffg=gff[which(gff$feature=="gene"),]
 
 out=data.frame()
 le=vector()
 li=vector()
 lASorf=vector()
 lASas=vector()

 #done=data.frame()
for (i in 1:nrow(gffg))
#for (i in 1:100)
 {

 repetitive="U"
#print(i)
#print(gffg$Name[i])

 if(gffg$Name[i]=="LTR")
 {
  name=paste("LTR",gffg$start[i],gffg$chr[i],sep='.')
 }
 
 else
 {
  name=gffg$Name[i]
 }

#print(name)

 ne=paste(name,"E",sep='.')
 ni=paste(name,"I",sep='.')
 ce=grep(ne, Chr[,1])
 te=grep(ne, Trs[,1])
 ci=grep(ni, Chr[,1])
 ti=grep(ni, Trs[,1])
 gc=which(gChr[,1]==name)
 gt=which(gTrs[,1]==name)

 re1=sum(Chr[ce, "reads"])
 re2=sum(Trs[te, "reads"])
 re=re1+re2
 ri1=sum(Chr[ci, "reads"])
 ri2=sum(Trs[ti, "reads"])
 ri=ri1+ri2
 le1=sum((Chr[ce,"FEAT.end"]-Chr[ce,"FEAT.start"]+1))/1000
 #le2=sum((Trs[te,"FEAT.end"]-Trs[te,"FEAT.start"]+1))/1000
 le[i]=le1
 li1=sum((Chr[ci,"FEAT.end"]-Chr[ci,"FEAT.start"]+1))/1000
 #li2=sum((Trs[ti,"FEAT.end"]-Trs[ti,"FEAT.start"]+1))/1000
 li[i]=li1

 if (!(name %in% gChr[,1]))
 {
  rAS=NA
  lASorf[i]=NA
  lASas[i]=NA
 }

 else
 {
  rAS1=gChr[gc,"AS.reads"]
  rAS2=gTrs[gt,"AS.reads"]
  rAS=rAS1+rAS2
  lASorf1=(gChr[gc,"ORF.end"]-gChr[gc,"ORF.start"]+1)/1000
  #lASorf2=(gTrs[gt,"ORF.end"]-gTrs[gt,"ORF.start"]+1)/1000
  lASorf[i]=lASorf1
    
  #print (gChr[which(gChr[,1]==name),])
  if (gChr[gc,"start"]=="NON_UNIQUE")
  {
   repetitive="M"
  }
  if (gChr[gc,"AS.start"]=="NON_UNIQUE")
  {
   repetitive="M"
   lASas1=NA
  }
  else
  {
   lASas1=(as.numeric(gChr[gc,"AS.end"])-as.numeric(gChr[gc,"AS.start"])+1)/1000
  }
#lASas2=(gTrs[which(gTrs[,1]==name),"AS.end"]-gTrs[which(gTrs[,1]==name),"AS.start"]+1)/1000
  lASas[i]=lASas1
  #if (le1 != le2){print("problem e...")}
  #if (li1 != li2){print("problem i...")}
  #if (lASorf1 != lASorf2){print("problem ASorf...")}
#if (lASas1 != lASas2){print("problem ASas...")}
 }
 
#cat(re,ri,rAS,"\n",sep='\t')
 out[i,1]=name
 out[i,2]=gffg$start[i]
 out[i,3]=gffg$end[i]
 out[i,4]=gffg$strand[i]
 out[i,5]=gffg$chr[i]
 out[i,6]=repetitive
 out[i,7]=re
 out[i,8]=ri
 out[i,9]=rAS
 out[i,10]=rAS
 out[i,11]=re
 out[i,12]=ri
 out[i,13]=rAS
 out[i,14]=0
 out[i,15]=0
 out[i,16]=0
 out[i,17]=0

 }
 
 #out[,7]=apply((out[,7]/le/depth),1,signif,digits=4)
 out[,7]=signif((out[,7]/le/depth), digits=4)
 #out[,8]=apply((out[,8]/li/depth),1,signif,digits=4)
 out[,8]=signif((out[,8]/li/depth), digits=4)
 #out[,9]=apply((out[,9]/lASorf/depth),1,signif,digits=4)
 out[,9]=signif((out[,9]/lASorf/depth), digits=4)
 #out[,10]=apply((out[,10]/lASas/depth),1,signif,digits=4)
 out[,10]=signif((out[,10]/lASas/depth), digits=4)
 out[,14]=le
 out[,15]=li
 out[,16]=lASorf
 out[,17]=lASas





colnames(out)=c("name","start","end","strand","chr","repetitive", "RPKM_exonic","RPKM_intronic","RPKM_AS","RPKM_ASexp","RC_exonic","RC_intronic","RC_AS","length_exonic", "length_intronic","length_orf","length_AS")
 return(out)
}




computeRPKM1<-function(f0, gff=gff, gff.name, depth)
{
 j=0
 Chr=read.delim(paste(f0,"col","mapIV",gff.name,"READSFEATURETABLE",sep='.'), stringsAsFactors=F)
 #Trs=read.delim(paste(f0,"ALLtrs","GEN","mapIV",gff.name,"READSFEATURETABLE",sep='.'), stringsAsFactors=F)
 gChr=read.delim(paste(f0,"col","mapIV",gff.name,"READSGENETABLE",sep='.'), stringsAsFactors=F)
 #gTrs=read.delim(paste(f0,"ALLtrs","GEN","mapIV",gff.name,"READSGENETABLE",sep='.'), stringsAsFactors=F)
 colnames(gChr)=c("Name","chr","strand","ORF.start","ORF.end","reads","start","end","AS.reads","AS.start","AS.end")
 #colnames(gTrs)=c("Name","chr","strand","ORF.start","ORF.end","reads","start","end","AS.reads","AS.start","AS.end")
 gffg=gff[which(gff$feature=="gene"),]
 
 out=data.frame()
 le=vector()
 li=vector()
 lASorf=vector()
 lASas=vector()

 #done=data.frame()
for (i in 1:nrow(gffg))
#for (i in 1:100)
 {

 repetitive="U"
#print(i)
#print(gffg$Name[i])

 if(gffg$Name[i]=="LTR")
 {
  name=paste("LTR",gffg$start[i],gffg$chr[i],sep='.')
 }
 
 else
 {
  name=gffg$Name[i]
 }

#print(name)

 ne=paste(name,"E",sep='.')
 ni=paste(name,"I",sep='.')
 ce=grep(ne, Chr[,1])
 #te=grep(ne, Trs[,1])
 ci=grep(ni, Chr[,1])
 #ti=grep(ni, Trs[,1])
 gc=which(gChr[,1]==name)
 #gt=which(gTrs[,1]==name)

 re1=sum(Chr[ce, "reads"])
 #re2=sum(Trs[te, "reads"])
 re=re1
 ri1=sum(Chr[ci, "reads"])
 #ri2=sum(Trs[ti, "reads"])
 ri=ri1
 le1=sum((Chr[ce,"FEAT.end"]-Chr[ce,"FEAT.start"]+1))/1000
 #le2=sum((Trs[te,"FEAT.end"]-Trs[te,"FEAT.start"]+1))/1000
 le[i]=le1
 li1=sum((Chr[ci,"FEAT.end"]-Chr[ci,"FEAT.start"]+1))/1000
 #li2=sum((Trs[ti,"FEAT.end"]-Trs[ti,"FEAT.start"]+1))/1000
 li[i]=li1

 if (!(name %in% gChr[,1]))
 {
  rAS=NA
  lASorf[i]=NA
  lASas[i]=NA
 }

 else
 {
  rAS1=gChr[gc,"AS.reads"]
  #rAS2=gTrs[gt,"AS.reads"]
  rAS=rAS1
  lASorf1=(gChr[gc,"ORF.end"]-gChr[gc,"ORF.start"]+1)/1000
  #lASorf2=(gTrs[gt,"ORF.end"]-gTrs[gt,"ORF.start"]+1)/1000
  lASorf[i]=lASorf1
    
  #print (gChr[which(gChr[,1]==name),])
  if (gChr[gc,"start"]=="NON_UNIQUE")
  {
   repetitive="M"
  }
  if (gChr[gc,"AS.start"]=="NON_UNIQUE")
  {
   repetitive="M"
   lASas1=NA
  }
  else
  {
   lASas1=(as.numeric(gChr[gc,"AS.end"])-as.numeric(gChr[gc,"AS.start"])+1)/1000
  }
#lASas2=(gTrs[which(gTrs[,1]==name),"AS.end"]-gTrs[which(gTrs[,1]==name),"AS.start"]+1)/1000
  lASas[i]=lASas1
  #if (le1 != le2){print("problem e...")}
  #if (li1 != li2){print("problem i...")}
  #if (lASorf1 != lASorf2){print("problem ASorf...")}
#if (lASas1 != lASas2){print("problem ASas...")}
 }
 
#cat(re,ri,rAS,"\n",sep='\t')
 out[i,1]=name
 out[i,2]=gffg$start[i]
 out[i,3]=gffg$end[i]
 out[i,4]=gffg$strand[i]
 out[i,5]=gffg$chr[i]
 out[i,6]=repetitive
 out[i,7]=re
 out[i,8]=ri
 out[i,9]=rAS
 out[i,10]=rAS
 out[i,11]=re
 out[i,12]=ri
 out[i,13]=rAS
 out[i,14]=0
 out[i,15]=0
 out[i,16]=0
 out[i,17]=0

 }
 
 #out[,7]=apply((out[,7]/le/depth),1,signif,digits=4)
 out[,7]=signif((out[,7]/le/depth), digits=4)
 #out[,8]=apply((out[,8]/li/depth),1,signif,digits=4)
 out[,8]=signif((out[,8]/li/depth), digits=4)
 #out[,9]=apply((out[,9]/lASorf/depth),1,signif,digits=4)
 out[,9]=signif((out[,9]/lASorf/depth), digits=4)
 #out[,10]=apply((out[,10]/lASas/depth),1,signif,digits=4)
 out[,10]=signif((out[,10]/lASas/depth), digits=4)
 out[,14]=le
 out[,15]=li
 out[,16]=lASorf
 out[,17]=lASas





colnames(out)=c("name","start","end","strand","chr","repetitive", "RPKM_exonic","RPKM_intronic","RPKM_AS","RPKM_ASexp","RC_exonic","RC_intronic","RC_AS","length_exonic", "length_intronic","length_orf","length_AS")
 return(out)
}




