anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/chicken_Gg6_atac_pre/"
tf=read.table("/scicore/home/tschoppp/GROUP/references/genomes/TFs/Gallus_gallus_TF_AnimalTFDB3_edited.txt", sep = "\t", header = T)
my.W = read.table("~/scRNA/ensembl97_wgenes.txt",sep="\t",header = TRUE,row.names = 1,stringsAsFactors = FALSE)

###################
## generate tfinfo
################### 

stamp2exp=function(subset.rna=subset.rna, idents="fine",
                   mf.each=mf.each, match.res=res,
                   rna.cls="sub_4",
                   atac.cls="sub_5"){
  
  #get average expression
  cells.rna=rownames(subset.rna@meta.data[subset.rna@meta.data[,idents] %in% rna.cls,])
  subset.cls=subset(subset.rna,cells=cells.rna)
  subset.cls$newset="Y"
  tmp=GetAssayData(subset.cls,assay = "RNA",slot="data")
  rownames(tmp)=gene.info[rownames(tmp),"genename"]
  subset.cls[['rename']]=CreateAssayObject(data=as.matrix(tmp))
  Idents(subset.cls)=subset.cls@meta.data[,"newset"]
  avg=AverageExpression(subset.cls,assays = "rename")[["rename"]]
  avg.cls=avg[,1]
  
  #get percentage exp
  #subset.cls=subset(subset.rna,idents = rna.cls)
  tmp=GetAssayData(subset.cls,assay = "RNA",slot = "counts")
  rownames(tmp)=gene.info[rownames(tmp),"genename"]
  pct.exp = rowSums(tmp>0)/ncol(tmp)*100
  
  #generate tf info matrix
  res1=match.res$pccCustom; #table(res1$searchMF)
  res1$motifNameNew=mid2name[res1$motifID,"motif_name"]
  n=6 #the number of predicted match for each de novo motif
  motifs=unique(res1$searchMF)
  tf.mtx2=matrix(nrow=length(motifs)*n,ncol=n)
  colnames(tf.mtx2)=c("motif.id","features.plot","features.name","avg.exp","pct.exp","negLog10Evalue")
  tf.mtx2[,"motif.id"]=rep(motifs,each=n)
  tf.mtx2[,"features.plot"]=rep(c(paste0("predict.tf",1:5),"predict.homerTF"),n=length(motifs))
  tf.mtx2=as.data.frame(tf.mtx2)
  
  #fill in tf info matrix
  for(i in 1:length(motifs)){
    tmp1=res1[res1$searchMF==motifs[i],"motifNameNew"]
    for(j in 1:length(tmp1)){
      tf.mtx2[tf.mtx2$motif.id==motifs[i] & tf.mtx2$features.plot==paste0("predict.tf",j),"features.name"]=toupper(tmp1[j])
      tf.mtx2[tf.mtx2$motif.id==motifs[i] & tf.mtx2$features.plot==paste0("predict.tf",j),"negLog10Evalue"]=res1[res1$searchMF==motifs[i],"negLog10_Evalue"][j]
    }
    
    tmp2.raw=res1[res1$searchMF==motifs[i],"bestMatchHomer"][1]
    if(!is.na(unlist(strsplit(tmp2.raw,"\\/"))[2])){
      tmp2.mod=paste0(unlist(strsplit(tmp2.raw,"\\/"))[1],"/",unlist(strsplit(tmp2.raw,"\\/"))[2])
    }else{
      tmp2.mod=paste0(unlist(strsplit(tmp2.raw,"\\/"))[1])
    }
    
    if(tmp2.mod %in% mid2name$motif_id2){
      tmp2.res=mid2name[mid2name$motif_id2 %in% tmp2.mod,"motif_name"]
    }else if(grepl("Jaspar",unlist(strsplit(tmp2.raw,"\\/"))[2]) | grepl("Jaspar",unlist(strsplit(tmp2.raw,"\\/"))[3])){
      if(grepl("MA",tmp2.mod)){
        tmp2.res=unlist(strsplit(tmp2.mod,"\\/"))[1]
      }else{
        tmp2.res=unlist(strsplit(tmp2.mod,"\\/|_"))[2]
      }
    }else{
      tmp2.res=mid2name[grepl(unlist(strsplit(tmp2.mod,"\\/|\\("))[1],mid2name$motif_id2),"motif_name"][1]
    }
    if(is.na(tmp2.res)){
      tmp2.res=unlist(strsplit(tmp2.mod,"\\/|\\("))[1]
    }
    tf.mtx2[tf.mtx2$motif.id==motifs[i] & tf.mtx2$features.plot=="predict.homerTF","features.name"]=toupper(tmp2.res)
    tf.mtx2[tf.mtx2$motif.id==motifs[i] & tf.mtx2$features.plot=="predict.homerTF","negLog10Evalue"]=NA
  }
  tf.mtx2$features.name=gsub("\\(VAR.2\\)|\\(VAR.3\\)","",tf.mtx2$features.name)
  tf.mtx2$negLog10Evalue=as.numeric(tf.mtx2$negLog10Evalue)
  
  tf.mtx2$avg.exp=avg.cls[tf.mtx2$features.name]
  tf.mtx2$avg.exp.log=log1p(tf.mtx2$avg.exp)
  tf.mtx2$pct.exp=pct.exp[tf.mtx2$features.name]
  
  #also consider the motifs mapped to several gene
  for(i in 1:nrow(tf.mtx2)){
    tmp=tf.mtx2[i,]
    if(grepl("\\|",tmp["features.name"])){
      tmp.gene=gsub("\\|",",",tmp["features.name"])
      tmp.gene=unlist(strsplit(tmp.gene,","))
      tmp.gene=tmp.gene[tmp.gene %in% names(avg.cls)]
      tf.mtx2[i,"avg.exp"]=mean(avg.cls[tmp.gene])
      tf.mtx2[i,"avg.exp.log"]=log1p(mean(avg.cls[tmp.gene]))
      tf.mtx2[i,"pct.exp"]=mean(pct.exp[tmp.gene])
    }
  }
  
  #get Homer pval
  if(length(atac.cls) > 1){
    mf.tmp=mf.each
  }else{
    mf.tmp=mf.each[[atac.cls]]
  }
  mf.info=matrix(ncol=3)
  for(i in 1:length(mf.tmp)){
    tmp.nm=names(mf.tmp)[i]
    tmp.mf=mf.tmp[[tmp.nm]]
    tmp.logP=tmp.mf@tags$logPval
    tmp.pval=unlist(strsplit(tmp.mf@tags$occurInfo,","))[3]
    tmp.pval=gsub("P:","",tmp.pval)
    mf.info=rbind(mf.info,c(tmp.nm,tmp.logP,tmp.pval))
  }
  mf.info=as.data.frame(mf.info[2:nrow(mf.info),]);dim(mf.info)
  colnames(mf.info)=c("motif.id","logP","pval")
  rownames(mf.info)=substr(mf.info$motif.id,1,20)
  tf.mtx2$homerPval=mf.info[tf.mtx2$motif.id,"pval"]
  
  #order by Evalue of predict.tf1
  tf.mtx2.topEvalue=tf.mtx2[tf.mtx2$features.plot=="predict.tf1",]
  tf.mtx2.topEvalue=tf.mtx2.topEvalue[order(-tf.mtx2.topEvalue$negLog10Evalue),]
  tf.mtx2$motif.id=factor(tf.mtx2$motif.id,levels = tf.mtx2.topEvalue$motif.id)
  
  return(tf.mtx2)
}


###################
## plot of TF expression & motif E-value
################### 

evalue.plot=function(tf.mtx=tf.mtx2){
  tf.mtx.topEvalue=tf.mtx[tf.mtx$features.plot=="predict.tf1",]
  tf.mtx.topEvalue$motif.id=factor(tf.mtx.topEvalue$motif.id, levels=levels(tf.mtx$motif.id))
  
  p1=ggplot(tf.mtx.topEvalue,aes(x=motif.id,y=negLog10Evalue))+theme_classic() +
    geom_point(col="red") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30))+ggtitle("top Evalue")
  
  
  p2=ggplot(tf.mtx,aes(x=motif.id,y=negLog10Evalue))+theme_classic() +
    geom_point(col="grey") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30)) +ggtitle("mean +/- sd")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")
  
  print(p1 + geom_hline(yintercept=seq(2,16,by=2),linetype="dotted"))
  print(p2 + geom_hline(yintercept=seq(2,16,by=2),linetype="dotted"))
}

rank.plot=function(tf.mtx=tf.cor$tfinfo,check="mf2exp.cor",
                   title="",
                   dotline.seq=seq(-0.5,0.5,by=0.25)){
  mfs=unique(tf.mtx$motif.id)
  top.mtx=NULL #matrix(nrow=length(mfs),ncol=ncol(tf.mtx))
  
  for(i in 1:length(mfs)){
    tmp=tf.mtx[tf.mtx$motif.id==mfs[i],]
    rank=tmp[!is.na(tmp[,check]),check]
    if(length(rank)==0){ 
      top.mtx=rbind(top.mtx,tmp[1,])  
    }else{
      top.mtx=rbind(top.mtx,tmp[!is.na(tmp[,check]) & tmp[,check]==max(rank),][1,])
    }
  }
  
  top.mtx=top.mtx[order(-top.mtx[,check]),]
  top.mtx$motif.id=factor(top.mtx$motif.id, levels=top.mtx$motif.id)
  
  top.plot=data.frame(id=top.mtx$motif.id,
                      value=top.mtx[,check])
  p1=ggplot(top.plot,aes(x=id,y=value))+ylab(check)+theme_classic() +
    geom_point(col="red") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30))+ggtitle(paste0(title,":top ",check," value"))
  
  tf.mtx$motif.id=factor(tf.mtx$motif.id, levels=top.mtx$motif.id)
  tf.plot=data.frame(id=tf.mtx$motif.id,
                     value=tf.mtx[,check])
  p2=ggplot(tf.plot,aes(x=id,y=value))+theme_classic() +
    geom_point(col="grey") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30)) +ggtitle("mean +/- sd")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")
  
  print(p1 + geom_hline(yintercept=dotline.seq,linetype="dotted"))
  print(p2 + geom_hline(yintercept=dotline.seq,linetype="dotted"))
}

cor.colors=c("#A50B89","#1207A3","#23AE22","#FFFEFE","#F0E420","#E63615","#682C36")
pval.colors=c("gray60","gray90","#FEF001","#FFCE03", "#FD9A01", "#FD6104", "#FF2C05", "#F00505")

tfmtx.plot=function(data.plot=tf.mtx2,cols.use = c("lightgrey", "blue"),dot.scale = 6,title="",
                    colfor="negLog10Evalue",colfor.use=c("blue","white","red"),
                    setbreaks=F,breaks=c(0,1),limits = c(-0.6,0.6)){
  plot1 <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'motif.id')) +
    geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradientn(colours=cols.use) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(x = 'Features', y = 'Identity') +
    theme_classic()
  
  if(setbreaks){
    plot2 <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'motif.id')) +
      geom_tile(aes_string(fill = colfor)) + 
      geom_text(aes(label = features.name),size=3) +
      #scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
      #scale_fill_gradientn(colours = c("gray90","gray90","gray80","yellow", "orange", "red", "darkred", "darkred")) +
      scale_fill_gradientn(colours = colfor.use,
                           values=scales::rescale(breaks),
                           limits =limits,breaks=breaks,na.value = "gray60") +
      theme_classic()
  }else{
    plot2 <- ggplot(data = data.plot, mapping = aes_string(x = 'features.plot', y = 'motif.id')) +
      geom_tile(aes_string(fill = colfor)) + 
      geom_text(aes(label = features.name),size=3) +
      scale_fill_gradientn(colours = colfor.use) +
      theme_classic()
  }
  
  print(plot1+ggtitle(title)+plot2)
}

###################
## z-score avgXpct
################### 

ztrans=function(tf.mtx=tf.mtx2){
  tf.mtx$avgXpct=tf.mtx$avg.exp*tf.mtx$pct.exp
  
  plot1 <- ggplot(data = tf.mtx, mapping = aes_string(x = 'features.plot', y = 'motif.id')) +
    geom_tile(aes_string(fill = 'avgXpct')) + 
    #geom_text(aes(label = features.name),size=3) +
    scale_fill_gradientn(colours=rainbow(7)[6:1]) + theme_classic()
  
  #center & scale the value
  tf.mtx.scale=matrix(ncol=6)
  colnames(tf.mtx.scale)=unique(tf.mtx$features.plot)
  for(i in 1:length(unique(tf.mtx$motif.id))){
    tmp.id=unique(tf.mtx$motif.id)[i]
    tmp.mtx=tf.mtx[tf.mtx$motif.id==tmp.id,]
    rownames(tmp.mtx)=tmp.mtx$features.plot
    tf.mtx.scale=rbind(tf.mtx.scale,tmp.mtx[colnames(tf.mtx.scale),"avgXpct"])
  }
  tf.mtx.scale=tf.mtx.scale[2:nrow(tf.mtx.scale),]
  rownames(tf.mtx.scale)=unique(tf.mtx$motif.id)
  tf.mtx.scale=t(scale(t(tf.mtx.scale)))
  tf.mtx.scale[is.nan(tf.mtx.scale)]=0
  
  tf.mtx.zmelt=melt(tf.mtx.scale)
  tf.mtx.zmelt$Var1=factor(tf.mtx.zmelt$Var1,levels = levels(tf.mtx$motif.id))
  tf.mtx.zmelt$Var2=factor(tf.mtx.zmelt$Var2,levels = c("predict.homerTF",unique(tf.mtx$features.plot)[1:5]))
  colnames(tf.mtx.zmelt)=c("motif.id","features.plot","zscore")
  plot2 <- ggplot(data = tf.mtx.zmelt, mapping = aes_string(x = 'features.plot', y = 'motif.id')) +
    geom_tile(aes_string(fill = 'zscore')) + 
    #geom_text(aes(label = features.name),size=3) +
    scale_fill_gradientn(colours=rainbow(7)[6:1]) + theme_classic()
  
  print(plot1+plot2)
  return(tf.mtx.scale)
}
  
ztrans.distriplot=function(tf.mtx=tf.mtx2,tf.mtx.scale=tf.mtx.scale){
  tf.mtx.scalemax=data.frame(motif.id=rownames(tf.mtx.scale),
                              zscoremax=apply(tf.mtx.scale,1,function(x){max(x[!is.na(x)])}))
  tf.mtx.scalemax=tf.mtx.scalemax[order(-tf.mtx.scalemax$zscoremax),]
  tf.mtx.scalemax$motif.id=factor(tf.mtx.scalemax$motif.id,levels = tf.mtx.scalemax$motif.id)
  p1=ggplot(tf.mtx.scalemax,aes(x=motif.id,y=zscoremax))+theme_classic() +
    geom_point(col="red") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30))+ggtitle("top Zscore")
  
  tf.mtx.zmelt=melt(tf.mtx.scale)
  tf.mtx.zmelt$Var1=factor(tf.mtx.zmelt$Var1,levels = levels(tf.mtx$motif.id))
  tf.mtx.zmelt$Var2=factor(tf.mtx.zmelt$Var2,levels = c("predict.homerTF",unique(tf.mtx$features.plot)[1:5]))
  colnames(tf.mtx.zmelt)=c("motif.id","features.plot","zscore")
  tf.mtx.zmelt$motif.id=factor(tf.mtx.zmelt$motif.id,levels = tf.mtx.scalemax$motif.id)
  tf.mtx.zmelt=tf.mtx.zmelt[order(tf.mtx.zmelt$motif.id),]
  p2=ggplot(tf.mtx.zmelt,aes(x=motif.id,y=zscore,na.rm=T))+theme_classic() +
    geom_point(col="grey") + theme(axis.text.x=element_text(angle=90, hjust=1,size = 10),axis.text.y=element_text(size = 25),axis.title.x=element_text(size = 30),axis.title.y=element_text(size = 30)) +ggtitle("mean +/- sd")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")
  
  print(p1 + geom_hline(yintercept=seq(0,max(tf.mtx.scalemax$zscoremax),by=0.5),linetype="dotted"))
  print(p2 + geom_hline(yintercept=seq(-max(tf.mtx.scalemax$zscoremax),max(tf.mtx.scalemax$zscoremax),by=0.5),linetype="dotted"))
}

###################
## generate aggregates/metacell
################### 
library(FNN)

embedding.perstage=function(rna=subset.rna,atac=subset.atac,usedrna = 'RNA',usedatac = 'activity',stages=c("L21","L24")){
  #run CCA with HVGs to find the nearest neighbor pairs::Matching of single-cell transcriptomes and epigenomes
  cpair=list()
  for(i in stages){
    srna=subset(rna,subset = orig.ident == i)
    satac=subset(atac,subset = dataset == i)
    srna=FindVariableFeatures(srna,assay=usedrna)
    hvgs=VariableFeatures(srna,assay=usedrna)[!(VariableFeatures(srna,assay=usedrna) %in% my.W$ensembl_gene_id)]
    hvgs = hvgs[hvgs %in% rownames(satac[[usedatac]])]
    print(length(hvgs))
    
    object.pair=NULL
    num.cc = 30
    if(ncol(srna) < 30 | ncol(satac) < 30){
      num.cc=min(ncol(srna),ncol(satac))-1
      print("Warning! the cell number of data is less than 30!")
    }
    object.pair = RunCCA(
      object1 = srna, object2 = satac,
      assay1 = usedrna, assay2 = usedatac,
      features = hvgs,
      num.cc = num.cc
    )
    object.pair <- L2Dim(object = object.pair, reduction = "cca")
    object.pair@meta.data$datatype = "RNA"
    object.pair@meta.data[rownames(satac@meta.data),"datatype"]="ATAC"
    p0=DimPlot(object = object.pair, reduction = "cca.l2", group.by = "datatype")+ggtitle(paste0(i,": CCA of scRNA & scATAC"))
    projs=Embeddings(object.pair,reduction="cca.l2")
    
    print(p0)
    
    cpair[[paste0(i)]]=projs
  }
  return(cpair)
}

library(ArchR)
#function
# https://github.com/GreenleafLab/ArchR/blob/968e4421ce7187a8ac7ea1cf6077412126876d5f/R/IntegrativeAnalysis.R
.computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
){
  .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  .requirePackage("nabor", source = "cran")
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}

determineOverlapCpp <- function(m, overlapCut) {
  .Call('_ArchR_determineOverlapCpp', PACKAGE = 'ArchR', m, overlapCut)
}

metapair=function(rna=mtx.fl_rna,atac=hm.integrated,
                  embeds=embeds,n=500, k=100){
  # 1) random sampling 500 scATAC cells as seed
  set.seed(1234)
  rD = Embeddings(atac,reduction = "lsi_bin")[,2:30]
  if(n > ncol(atac)){
    seed=sample(rownames(atac@meta.data),n,replace = T)
  }else{
    seed=sample(rownames(atac@meta.data),n)
  }
  knnObj <- knnx.index(data=rD, query=rD[seed,], k=k,algorithm = "kd_tree");dim(knnObj)
  # 2) for each seed, get its 99 nearest scATAC neighbors, and remove if the one with > 80% cells overlapped
  overlapCutoff=0.8
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  knnObj <- knnObj[keepKnn==0,];dim(knnObj)
  if(nrow(knnObj) < 50){
    print(paste0("Warning! only get few low-overlapping aggregates!",nrow(knnObj)))
  }
  # 3) then in CCA space find the 100 nearest scRNA neighbors of seed (use the result from previous CCA)
  projs=embeds
  cell.rna=rownames(rna@meta.data)
  index=knnx.index(data=projs[cell.rna,], query=projs[seed[keepKnn==0],], k=k,algorithm = "kd_tree");dim(index)
  aggr=list()
  aggr[["atac.cname"]]=rownames(rD); aggr[["atac.corder"]]=knnObj
  aggr[["rna.cname"]]=cell.rna; aggr[["rna.corder"]]=index
  return(aggr)
  
}

metapair.perstage=function(allrna=subset.rna,allatac=subset.atac,
                           embeds=embeds,n=100, k=20,stages=c("L21","L24")){
  # 1) random sampling scATAC cells as seed
  aggr=list()
  for(i in stages){
    rna=subset(allrna,subset = orig.ident == i)
    atac=subset(allatac,subset = dataset == i)
    set.seed(1234)
    rD = Embeddings(atac,reduction = "lsi_bin")[,2:30]
    seed=sample(rownames(atac@meta.data),n)
    knnObj <- knnx.index(data=rD, query=rD[seed,], k=k,algorithm = "kd_tree");dim(knnObj)
    # 2) for each seed, get its 99 nearest scATAC neighbors, and remove if the one with > 80% cells overlapped
    overlapCutoff=0.8
    keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
    knnObj <- knnObj[keepKnn==0,];print(dim(knnObj))
    # 3) then in CCA space find the 100 nearest scRNA neighbors of seed (use the result from previous CCA)
    projs=embeds[[paste0("embedding_",i)]]
    cell.rna=rownames(rna@meta.data)
    index=knnx.index(data=projs[cell.rna,], query=projs[seed[keepKnn==0],], k=k,algorithm = "kd_tree");dim(index)
    
    aggr[[i]][["atac.cname"]]=rownames(rD); aggr[[i]][["atac.corder"]]=knnObj
    aggr[[i]][["rna.cname"]]=cell.rna; aggr[[i]][["rna.corder"]]=index
    
  }
  return(aggr)
}


###################
## calculate motif activity
################### 
library(BSgenome.Ggallus.ensembl.galGal6)

toPFMlist=function(mf.ppmtx=mf.tmp){
  mf2pfm=list()
  for(i in 1:length(mf.ppmtx)){
    tmp.mtx=PFMatrix(ID=mf.ppmtx[[i]]@ID, name=mf.ppmtx[[i]]@name, 
                     bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                     tags=list(species="Gallus gallus",
                               logOdds=mf.ppmtx[[i]]@tags$logOdds, logPval=mf.ppmtx[[i]]@tags$logPval,
                               occurInfo=mf.ppmtx[[i]]@tags$occurInfo,
                               tax_group="vertebrates"),
                     profileMatrix=round(mf.ppmtx[[i]]@profileMatrix*20,0)
    )
    mf2pfm[[mf.ppmtx[[i]]@ID]]=tmp.mtx
  }
  mfList=do.call(PFMatrixList,mf2pfm)
  return(mfList)
}

###################
## calculate sum activity & expression
################### 

metamtx=function(rna=subset.rna,atac=subset.atac,slot.rna="scale.data",
                 rna.assay="RNA",atac.assay="chromvar",stage=stages,
                 aggr=metapairs,genelist=unique(tf.mtx2$features.name)){
  
  #expression
  exp=GetAssayData(rna,assay = rna.assay,slot = slot.rna)
  exp=exp[rownames(exp) %in% rownames(gene.info),]
  rownames(exp)=gene.info[rownames(exp),"genename"]
  #also consider the motifs mapped to several gene
  for(i in 1:length(genelist)){
    tmp=genelist[i]
    if(grepl("\\|",tmp)){
      tmp.gene=gsub("\\|",",",tmp)
      tmp.gene=unlist(strsplit(tmp.gene,","))
      tmp.gene=tmp.gene[tmp.gene %in% rownames(exp)]
      if(length(tmp.gene)>1){
        exp=rbind(exp,apply(exp[tmp.gene,],2,mean))
        rownames(exp)[nrow(exp)]=tmp
      }else if(length(tmp.gene)==1){
        exp=rbind(exp,exp[tmp.gene,])
        rownames(exp)[nrow(exp)]=tmp
      }
    }
  }
  #only keep the ones in genelist
  genes.keep=rownames(exp)
  genes.keep=genes.keep[genes.keep %in% genelist]
  exp=exp[genes.keep,]
  aggr.expmtx=matrix(nrow=nrow(exp),ncol=nrow(aggr$atac.corder))
  for(i in 1:ncol(aggr.expmtx)){
    cells=aggr[["rna.cname"]][aggr[["rna.corder"]][i,]]
    aggr.expmtx[,i]=colSums(t(exp[,cells]))
  }
  rownames(aggr.expmtx)=rownames(exp)
  colnames(aggr.expmtx)=paste0(stage,"-","metacell",1:nrow(aggr$atac.corder))
  
  #motif activity deviation
  act=GetAssayData(atac,assay=atac.assay,slot = "data")
  aggr.actmtx=matrix(nrow=nrow(act),ncol=nrow(aggr$atac.corder))
  for(i in 1:ncol(aggr.actmtx)){
    cells=aggr[["atac.cname"]][aggr[["atac.corder"]][i,]]
    aggr.actmtx[,i]=colSums(t(act[,cells]))
  }
  rownames(aggr.actmtx)=substr(rownames(act),1,20)
  colnames(aggr.actmtx)=paste0(stage,"-","metacell",1:nrow(aggr$atac.corder))
  
  aggr.mtx=list()
  aggr.mtx[["expression"]]=aggr.expmtx
  aggr.mtx[["motifactivity"]]=aggr.actmtx
  return(aggr.mtx)
}

###################
## whole workflow
################### 

calMotifsCor=function(data.atac=hm.integrated,
                      data.rna=mtx.fl_rna,slot.rna="scale.data",
                      cca.rna="RNA", cca.atac="activity",
                      ident.rna="annotation_broad",
                      ident.atac="annotation_broad",
                      ctype.rna=c("Mesenchyme","skeletogenicCells"),
                      ctype.atac=c("Mesenchyme","skeletogenicCells","MesenchymeORskeletogenicCells"),
                      motifs=mf.each$sub_5, runchromvar=F,
                      mtx.rna="RNA", mtx.atac="chromvar",
                      stages=c("L21","L24"),
                      tfinfo=tf.mtx2,
                      cor.method="pearson",
                      k=20,n=200,hvgsN=3000
){
  aggr.mtx=list()
  for(x in 1:length(stages)){
    ## subset RNA & ATAC
    Idents(data.rna)=data.rna@meta.data[,ident.rna]
    cells.rna=rownames(data.rna@meta.data[data.rna@meta.data[,ident.rna] %in% ctype.rna & data.rna@meta.data[,"orig.ident"] %in% stages[x],])
    subset.rna=subset(data.rna,cells=cells.rna)
    
    Idents(data.atac)=data.atac@meta.data[,ident.atac]
    cells.atac=rownames(data.atac@meta.data[data.atac@meta.data[,ident.atac] %in% ctype.atac & data.atac@meta.data[,"dataset"] %in% stages[x],])
    subset.atac=subset(data.atac,cells=cells.atac)
    print("finish subsetting data!")
    
    ## aggregate pairs
    embeds=embedding.perstage(rna=subset.rna,atac=subset.atac,usedrna = cca.rna,usedatac = cca.atac,stages=stages[x])
    metapairs=metapair(rna=subset.rna,atac=subset.atac,
                       embeds=embeds[[stages[x]]],n=n, k=k)
    print("finish get aggregates!")
    
    ## calculate motif activity
    if(runchromvar){
      mflist=toPFMlist(mf.ppmtx=motifs)
      DefaultAssay(subset.atac)="macs2peaks"
      subset.atac <- AddMotifs(
        object = subset.atac, pfm = mflist,
        genome = BSgenome.Ggallus.ensembl.galGal6
      )
      subset.atac=RunChromVAR(subset.atac,genome=BSgenome.Ggallus.ensembl.galGal6) 
    }
    print("finish calculating motif activity!")
    
    ## restrict expression to highly variable ones
    subset.rna@active.assay=mtx.rna
    subset.rna=FindVariableFeatures(subset.rna,assay=mtx.rna,nfeatures=hvgsN)
    subset.rna=ScaleData(subset.rna,assay=mtx.rna,
                         features = VariableFeatures(subset.rna,assay=mtx.rna))
    #no need for motif activity, as we only care the dozens of de novo in specific cell type
    
    ## calculate sum of motif activity & expression
    aggr.mtx[[stages[x]]]=metamtx(rna=subset.rna,atac=subset.atac,slot.rna="scale.data",
                                  rna.assay=mtx.rna,atac.assay=mtx.atac,stage=stages[x],
                                  aggr=metapairs,genelist=unique(tfinfo$features.name))
    print("finish calculating sum & will run the next sample!")
    
  }
  print("finish running all samples!")
  
  ## put all aggregates data together
  aggr.all=list()
  aggr.all[["expression"]]=NULL;aggr.all[["motifactivity"]]=NULL
  names.act=rownames(aggr.mtx[[stages[1]]]$motifactivity)
  names.exp=NULL
  for(i in 1:length(stages)){
    names.exp=c(names.exp,rownames(aggr.mtx[[stages[i]]]$expression))
  }
  names.exp=names.exp[!duplicated(names.exp)]
  for(i in stages){
    #need to deal with non-overlap genes, put NA
    tmp=NULL
    tmp=aggr.mtx[[i]][["expression"]]
    tmp.lost=names.exp[!(names.exp %in% rownames(tmp))]
    tmp.add=matrix(nrow=length(tmp.lost),ncol=ncol(tmp))
    rownames(tmp.add)=tmp.lost;colnames(tmp.add)=colnames(tmp)
    tmp.new=rbind(tmp,tmp.add)
    
    aggr.all[["expression"]]=cbind(aggr.all[["expression"]],tmp.new[names.exp,])
    aggr.all[["motifactivity"]]=cbind(aggr.all[["motifactivity"]],aggr.mtx[[i]][["motifactivity"]][names.act,])
  }
  print("finally, start to calculate the correlation!")
  
  ## calculate correlation
  tfinfo$mf2exp.cor=NA
  tfinfo$mf2exp.corpval=NA
  tfinfo$motif.id=gsub("_","-",tfinfo$motif.id)
  for(i in 1:nrow(tfinfo)){
    if(tfinfo$features.name[i] %in% rownames(aggr.all$expression)){
      mfact.tmp=aggr.all$motifactivity[tfinfo$motif.id[i],]
      exp.tmp=aggr.all$expression[tfinfo$features.name[i],]
      #need to deal with  NA
      #will be omited by default
      cor=cor.test(x=exp.tmp,y=mfact.tmp,method = cor.method)
      tfinfo$mf2exp.cor[i]=cor$estimate
      tfinfo$mf2exp.corpval[i]=cor$p.value
    }
  }
  
  res=list()
  res[["aggr"]]=aggr.all
  res[["tfinfo"]]=tfinfo
  res[["hvgs"]]=names.exp
  return(res)
  
}

onlyCor=function(tfinfo=tf.mtx2,aggr.all=aggr.all,
                 cor.method="pearson"){
  ## calculate correlation
  tfinfo$mf2exp.cor=NA
  tfinfo$mf2exp.corpval=NA
  tfinfo$motif.id=gsub("_","-",tfinfo$motif.id)
  for(i in 1:nrow(tfinfo)){
    if(tfinfo$features.name[i] %in% rownames(aggr.all$expression)){
      mfact.tmp=aggr.all$motifactivity[tfinfo$motif.id[i],]
      exp.tmp=aggr.all$expression[tfinfo$features.name[i],]
      cor=cor.test(x=exp.tmp,y=mfact.tmp,method = cor.method)
      tfinfo$mf2exp.cor[i]=cor$estimate
      tfinfo$mf2exp.corpval[i]=cor$p.value
    }
  }
  return(tfinfo)
}

biplot=function(tf.cor=tf.cor,
                gene="SOX9",motif="clssub-5-1-RAACAAAGV",
                mfactTo0.1=F){
  if(!(gene %in% rownames(tf.cor$aggr$expression))){ 
    print(paste0("this gene is not HVGs:", gene))
    cormtx=NULL
  }else{
    cormtx=data.frame(expression=(tf.cor$aggr$expression[gene,]),
                      activity=(tf.cor$aggr$motifactivity[motif,]),
                      aggr=colnames(tf.cor$aggr$expression))
    cormtx=cormtx[!is.na(cormtx$expression),]
    if(mfactTo0.1){
      cormtx$activity=cormtx$activity/10
    }
    print("quantile expression");print(quantile(cormtx$expression))
    print("quantile activity");print(quantile(cormtx$activity))
    
    coeff=tf.cor$tfinfo[tf.cor$tfinfo$motif.id==motif & tf.cor$tfinfo$features.name==gene,"mf2exp.cor"]
    pval=tf.cor$tfinfo[tf.cor$tfinfo$motif.id==motif & tf.cor$tfinfo$features.name==gene,"mf2exp.corpval"]
    
    p1=ggplot(cormtx,aes(x=activity,y=expression, na.rm = TRUE)) + 
      ggtitle(paste0(motif," & ",gene,", coeff: ",round(coeff,2),", -log10(pval):",round(-log10(pval),2))) +
      geom_point()+ #aes(colour = species)
      xlim(min(cormtx$activity)-5,max(cormtx$activity)+5)+
      ylim(min(cormtx$expression)-5,max(cormtx$expression)+5)+
      theme_classic()+
      #geom_abline(slope=1, linetype="dashed",color="red")+
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
      xlab("motif activity")+ylab("gene expression") 
    print(p1)
  }
  
  return(cormtx)
}

###################
## assignment of motif name
################### 

# to include homer, preset negLog10Evalue=3

# for those with more than one candidate left:
# if there is significant mf2exp.cor, then keep the one with highest correlation coefficient
  # if not, all from STAMP, then keep the one with lowest E value
    # if not, keep the lowest E-value from STAMP, compare with the one from homer, keep the one with highest pct.exp

# if several motifs give same TF name, then keep only the one with smallest homer Pval

mf2names=function(tfscore=tf.cor$tfinfo, compareHomer="pct.exp",
                  negLogEval.cutoff=2, avgExp.cutoff=0, pctExp.cutoff=3,
                  rmDup=F){
  #preset negLog10Evalue=3
  tfscore=tfscore[!is.na(tfscore$motif.id),]
  tfscore[tfscore$features.plot=="predict.homerTF","negLog10Evalue"]=(negLogEval.cutoff+1)
  #general cutoff
  tfscore=tfscore[tfscore$negLog10Evalue > negLogEval.cutoff & (tfscore$avg.exp > avgExp.cutoff & tfscore$pct.exp > pctExp.cutoff),]
  
  mfid=unique(tfscore$motif.id);mfid=mfid[!is.na(mfid)]
  mfname=data.frame(motifid=mfid,
                    tfname=rep("unknown",length(mfid)),
                    homerPval=rep(1,length(mfid)),
                    defineClass=rep("classN",length(mfid)),
                    trend=rep("unknown",length(mfid)))
  #filter step
  for(i in 1:nrow(mfname)){
    tmp=tfscore[tfscore$motif.id==mfid[i],]
    tmp=tmp[!is.na(tmp$motif.id),]
    tmp[is.na(tmp$mf2exp.cor),"mf2exp.cor"]=0
    tmp[is.na(tmp$mf2exp.corpval),"mf2exp.corpval"]=1
    mfname$homerPval[i]=tmp$homerPval[1]
    tmp.sig=tmp[tmp$mf2exp.corpval<0.05,] # & tmp$mf2exp.cor>0
    if(nrow(tmp.sig)>0){
      mfname$tfname[i]=tmp.sig[tmp.sig$mf2exp.corpval==min(tmp.sig$mf2exp.corpval),"features.name"][1]
      mfname$defineClass[i]="class1"
      tmp.trend=tmp.sig[tmp.sig$mf2exp.corpval==min(tmp.sig$mf2exp.corpval),"mf2exp.cor"][1]
      if(tmp.trend>0){
        mfname$trend[i]="enhancer"
      }else{
        mfname$trend[i]="repressor"
      }
    }else if(!("predict.homerTF" %in% tmp$features.plot)){
      mfname$tfname[i]=tmp[tmp$negLog10Evalue==max(tmp$negLog10Evalue),"features.name"][1]
      mfname$defineClass[i]="class3"
    }else{
      keep=tmp[tmp$features.plot=="predict.homerTF" | tmp$negLog10Evalue==max(tmp$negLog10Evalue),]
      mfname$tfname[i]=keep[keep[,compareHomer]==max(keep[,compareHomer]),"features.name"][1]
      mfname$defineClass[i]="class2"
    }
    
  }
  print(paste0(nrow(mfname)," motifs are retained with names."))
  
  # if several motifs give same TF name, then keep only the one with smallest homer Pval
  if(rmDup){
    mfname=rmDuplicates(mfname=mfname)
  }
  
  return(mfname)
}

rmDuplicates=function(mfname=mf2tfnames){
  dup=unique(mfname$tfname[duplicated(mfname$tfname)])
  for(i in 1:length(dup)){
    tmp=mfname[mfname$tfname==dup[i],]
    tmp=tmp[order(tmp$homerPval),]
    mfname[mfname$motifid %in% tmp$motifid[2:nrow(tmp)],"tfname"]="duplicated"
  }
  mfname=mfname[mfname$tfname != "duplicated",]
  print(paste0("after remove motifs with same name, ",nrow(mfname)," left."))
  return(mfname)
}

save_pheatmap_pdf <- function(x, filename, width = 15,height = 12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

toSeqLogo=function(motifs=mf.each,id="sub_5.clssub_5_1_RAACAAAGVNNBCWTTGTT"){
  tmp.pwm=motifs[[id]]@profileMatrix
  tmp.pfm=as.matrix(round(tmp.pwm*20,0))
  rownames(tmp.pfm)=c("A","C","G","T")
  tmp.mtx=PFMatrix(ID=id,  
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   tags=list(species="Gallus gallus",
                             tax_group="vertebrates"),
                   profileMatrix=tmp.pfm
  )
  
  #seqLogo(toICM(tmp.mtx))
  return(tmp.mtx)
}

