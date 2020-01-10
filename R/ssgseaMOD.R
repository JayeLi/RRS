mod.generate <- function(data,md.nm='Default',md.dscrp='None',md.nt='rank',statistic='area.under.RES',weight=0.25,random.seed=12345,nperm=0,col=c('black','lightgrey'),link.func='logit',filter=c(0.5,-0.5),out='tmp.mod'){
  if(is.null(names(data))){
    tmp = data;
    data = rep(filter[1]+1,times=length(tmp))
    names(data) = tmp
  }
  up.genes = names(data)[data>filter[1]]
  dn.genes = names(data)[data<filter[2]]
  write.table(t(c("model.creation.date",paste("'",date(),"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F)
  write.table(t(c("model.name",paste("'",md.nm,"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("model.description",paste("'",md.dscrp,"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("sample.norm.type",paste("'",md.nt,"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("statistic",paste("'",statistic,"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("weight",weight)),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("random.seed",random.seed)),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("nperm",nperm)),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("link.function",paste("'",link.func,"'",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  write.table(t(c("cl",paste("c('",paste(col,collapse="','"),"')",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  if(length(up.genes)>0){
    write.table(t(c("msig.up.genes",paste("c('",paste(up.genes,collapse="','"),"')",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  }
  if(length(dn.genes)>0){
    write.table(t(c("msig.dn.genes",paste("c('",paste(dn.genes,collapse="','"),"')",sep=''))),out,sep='\t',quote=F,col.names=F,row.names=F,append=T)
  }
}


mod.analyze <- function(mat,mod.name,mod.dir='./',norm=T,all=F,out='tmp.gct'){
  library(ssgsea.GBM.classification)
  write.table('#1.2',out,quote=F,col.name=F,row.name=F)
  write.table(paste(nrow(mat),ncol(mat),sep="\t"),out,quote=F,col.name=F,row.name=F,append=T)
  write.table(t(c(c('NAME','DESCRIPTION'),colnames(mat))),out,quote=F,col.name=F,row.name=F,append=T,sep="\t")
  write.table(cbind(rownames(mat),rep('NA',times=nrow(mat)),mat),out,quote=F,col.name=F,row.name=F,append=T,sep='\t')
  OPAM.apply.model.2(
        input.ds           = out,
        models.dir         = mod.dir,
        models             = mod.name,
        raw.score.outfile  = paste(out,'.score',sep=''),
        norm.score.outfile = paste(out,'.norm.score',sep=''),
	model.score.outfile= paste(out,'.model.score',sep=''),
	prob.outfile= "",graphics.off= T)
  if(all){
  	mod.rs1 <- read.table(paste(out,'.norm.score',sep=''),skip=2,header=T,row.names=1,sep='\t')
  	mod.rs2 <- read.table(paste(out,'.score',sep=''),skip=2,header=T,row.names=1,sep='\t')
  	return( list(norm=mod.rs1[,-1],raw=mod.rs2[,-1]) )
  }else{
  	fm = ifelse(norm,'.norm.score','.score')
  	mod.rs <- read.table(paste(out,fm,sep=''),skip=2,header=T,row.names=1,sep='\t')
  	mod.rs <- mod.rs[,-1]
  	return(mod.rs)
  }
}


mod.analyze2 <- function(mat,mod.name,mod.dir='./',norm=F,all=T,permN=10000,out='tmp.gct'){

  random_profile<-c()
  for (i in 1:permN){
    a<-mat[cbind(seq(1:nrow(mat)),sample(1:ncol(mat),nrow(mat),replace=T))]  # the first column is the decription column, ignore it.
    random_profile<-cbind(random_profile,a)
    if(i%%100==0){
      print(i)
    }
  }
  rownames(random_profile) <- rownames(mat)

  random = t(mod.analyze(random_profile,mod.name,mod.dir,norm,all,out)$raw)
  data   = mod.analyze(mat,mod.name,mod.dir,norm,all,out)
  cat("[1] \"MOD Complete\"\n")
  tdata  = t(data$raw)
  p_result<-tdata;
  for (i in 1:dim(tdata)[1]){
    p_result[i,]<-colSums(sweep(random,2,tdata[i,])>=0);
  }
  cat(dim(p_result),'\n')
  colnames(p_result)=paste(colnames(p_result),"pval",sep="_");
  p_result=t(apply(p_result,1,function(x){(x+1)/(permN+1)}))
  cat('TWO\n')
  tdata = cbind(tdata,t(data$norm))
  cat("THREE\n")
  cat(dim(tdata),dim(p_result),"\n")
  colnames(tdata) = paste(rep(mod.name,times=2),rep(c("raw","norm"),times=c(length(mod.name),length(mod.name))),sep="_")
  if(length(mod.name)==1){
    tmp = cbind(tdata,t(p_result))
    colnames(tmp) = c(colnames(tdata),paste(mod.name,'pval',sep='_'))
    return(tmp)
  }
  return(cbind(tdata,p_result))

}
