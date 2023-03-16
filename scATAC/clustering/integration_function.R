## save own function

###################
## label transferring
###################

getpredict = function(atac=hm.integrated,rna=mtx.fl_rna,usedatac = 'activity', usedrna = 'RNA',dims.anchor=1:30,embed.transfer="harmony"){
  DefaultAssay(atac) <- usedatac; DefaultAssay(rna) <- usedrna 
  hvgs=VariableFeatures(rna)[!(VariableFeatures(rna) %in% my.W$ensembl_gene_id)]
  print(length(hvgs))
  hvgs = hvgs[hvgs %in% rownames(atac[[usedatac]])]
  print(length(hvgs))
  transfer.anchors <- FindTransferAnchors(
    reference = rna,query = atac,dims = dims.anchor,
    reference.assay = usedrna, query.assay = usedatac,
    max.features = 200, k.filter = 200, #default
    reduction = 'cca',features = hvgs
  )
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = Idents(rna),
    weight.reduction = atac[[embed.transfer]],
    dims = 2:30
  )
  return(predicted.labels)
}

stackedbar = function(atac=hm.integrated,input=prediction,thres=0.5,method="ArchR",rcolor=rcolor){
  predict = cbind(input[,c("predicted.id","prediction.score.max")],atac@meta.data[rownames(input),"cls"])
  colnames(predict) = c("predictedId","predictionScoreMax","cluster")
  predict[predict$predictionScoreMax <= thres,"predictedId"]="NA"; #str(predict)
  predict$predictedId=factor(predict$predictedId,levels = names(rcolor))
  g = ggplot(predict,aes(cluster)) + theme_classic()
  g1=g + geom_bar(aes(fill = predictedId)) +
    scale_fill_manual("legend", values = rcolor)
  g2=g + geom_bar(aes(fill = predictedId),position="fill") +
    geom_hline(yintercept=0.25, linetype="dashed", color = "red") +
    scale_fill_manual("legend", values = rcolor)
  print(g1+ggtitle(method)+g2+ggtitle("percentage"))
}

###################
## NNLS-non-negative least square
###################
#https://www.nature.com/articles/s41586-019-0969-x
#https://www.r-bloggers.com/2019/11/non-negative-least-squares/

getMtpNnls=function(select=daptop.all,test.rna=avg.rna,test.atac=avg.atac,norm=F){
  #scale
  test.rna=t(scale(t(test.rna)))
  test.atac=t(scale(t(test.atac)))
  #calculate beta, scRNA as target
  coeff.rna=NULL
  for(i in colnames(test.rna)){
    #nnls
    x=test.atac[select,]
    y=test.rna[select,i]
    #tmp.nnls=nnls(x,y)
    #tmp.coeff=tmp.nnls$x
    tmp.nnls=glmnet(x, y, lambda = 0, lower.limits = 0, intercept = FALSE)
    tmp.coeff=as.numeric(coef(tmp.nnls))[2:(ncol(x)+1)]
    coeff.rna=rbind(coeff.rna,tmp.coeff)
  }
  rownames(coeff.rna)=colnames(test.rna)
  colnames(coeff.rna)=colnames(test.atac)
  print(round(coeff.rna,2))
  coeff.rna[coeff.rna > 10]=0 #exclude extreme wierd value
  
  #calculate beta, scATAC as target
  coeff.atac=NULL
  for(i in colnames(test.atac)){
    #nnls
    x=test.rna[select,]
    y=test.atac[select,i]
    #tmp.nnls=nnls(x,y)
    #tmp.coeff=tmp.nnls$x
    tmp.nnls=glmnet(x, y, lambda = 0, lower.limits = 0, intercept = FALSE)
    tmp.coeff=as.numeric(coef(tmp.nnls))[2:(ncol(x)+1)]
    coeff.atac=rbind(coeff.atac,tmp.coeff)
  }
  rownames(coeff.atac)=colnames(test.atac)
  colnames(coeff.atac)=colnames(test.rna)
  print(round(coeff.atac,2))
  coeff.atac[coeff.atac > 10]=0 #exclude extreme wierd value
  
  #--> multiply the beta
  coeff=(coeff.rna+0.01)*(t(coeff.atac+0.01))
  
  #normalize coeff
  if(norm){
    for(i in 1:ncol(coeff)){
      coeff[,i]=((coeff[,i])/max(coeff[,i]))
    }
  }
  
  return(coeff)
}

###################
## summary barplot
###################

stackedbar.new = function(input=hm.integrated@meta.data,cname=c("dataset","annonew")){
  predict = input[,cname]
  colnames(predict) = c("stage","annotation")
  
  predict$stage=factor(predict$stage)
  g = ggplot(predict,aes(annotation)) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=1))
  g1=g + geom_bar(aes(fill = stage)) #+
  #scale_fill_manual("legend", values = rcolor)
  g2=g + geom_bar(aes(fill = stage),position="fill") #+
  #geom_hline(yintercept=0.25, linetype="dashed", color = "red") +
  #scale_fill_manual("legend", values = rcolor)
  
  perc=table(predict$stage,predict$annotation)  
  perc=perc/rowSums(perc)*100; mperc=melt(perc)
  tperc=t(perc);tperc=tperc/rowSums(tperc)*100; mperc=melt(t(tperc))
  g3=ggplot(mperc, aes(x = Var2, y = value)) + theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=1))+
    geom_col(aes(fill = Var1), width = 0.7) #+
  #scale_fill_manual("legend", values = cols)
  
  print(g1+ggtitle("cell numbers")+ #g2+ggtitle("percentage")+
          g3+ggtitle("percentage (normalize by stage)"))
  
}

###################
## xtSNE
###################
source("~/software/FIt-SNE/fast_tsne.R")
runxtsne=function(my.se=my.subset,my.dimensions=c(2:30),perplexity=30,perplexity_list = c(30, round(dim(my.subset)[2]/100)),origin.dims="harmony"){
  if(round(dim(my.se)[2]/12) < 200){
    learning_rate = 200
  }else{
    learning_rate = round(dim(my.se)[2]/12)
  }
  
  #xtsne, with exaggeration
  my.tsne = fftRtsne(my.se@reductions[[origin.dims]]@cell.embeddings[,my.dimensions],
                     max_iter = 1000,
                     perplexity = perplexity,
                     learning_rate =  learning_rate,
                     initialization =(my.se@reductions[[origin.dims]]@cell.embeddings[,1:2]/my.se@reductions[[origin.dims]]@stdev[1])*0.0001,
                     perplexity_list = perplexity_list,
                     late_exag_coeff = 4,
                     start_late_exag_iter = 250,
                     fast_tsne_path="~/software/FIt-SNE/bin/fast_tsne")
  rownames(my.tsne) = colnames(my.se)
  my.se[["xtsne"]] = CreateDimReducObject(embeddings = my.tsne, key = "xtsne_", assay = DefaultAssay(my.se), global = T)
  my.se = AddMetaData(my.se, metadata = data.frame(Embeddings(my.se, "xtsne")))
  
  #tsne
  my.tsne = fftRtsne(my.se@reductions[[origin.dims]]@cell.embeddings[,my.dimensions],
                     max_iter = 1000,
                     perplexity = perplexity,
                     learning_rate =  learning_rate,
                     initialization =(my.se@reductions[[origin.dims]]@cell.embeddings[,1:2]/my.se@reductions[[origin.dims]]@stdev[1])*0.0001,
                     perplexity_list = perplexity_list,
                     fast_tsne_path="~/software/FIt-SNE/bin/fast_tsne")
  rownames(my.tsne) = colnames(my.se)
  my.se[["tsne"]] = CreateDimReducObject(embeddings = my.tsne, key = "tsne_", assay = DefaultAssay(my.se), global = T)
  my.se = AddMetaData(my.se, metadata = data.frame(Embeddings(my.se, "tsne")))
  return(my.se)
}
