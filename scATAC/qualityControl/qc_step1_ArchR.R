library(ArchR); packageVersion("ArchR")
library(GenomeInfoDb); library(ensembldb)
set.seed(1234)

sample.name <- c("L21","L24","N15","N18","S12","S15")
data.name <- c("10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_12012021")
anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/chicken_Gg6_atac_pre/"

for(i in 1:length(sample.name)){
  data_path <- paste0("/scicore/home/tschoppp/GROUP/mapped_data/",data.name[i],"/",sample.name[i],"/outs/")
  
  ## Creating Arrow Files & ArchRProject
  ArrowFiles <- createArrowFiles(
    inputFiles = paste0(data_path,'fragments.tsv.gz'), 
    sampleNames = sample.name[i],
    geneAnnotation = readRDS(paste0(anno_path,"geneAnnotation.rds")),
    genomeAnnotation = readRDS(paste0(anno_path,"genomeAnnotation.rds")),
    minFrags = 100,
    filterTSS = 1, 
    filterFrags = 1000, 
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE
  )
  
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles , 
    outputDirectory = paste0("/scicore/home/tschoppp/wang0007/scATAC/reslt_scATAC/clustering/",sample.name[i]),
    geneAnnotation = readRDS(paste0(anno_path,"geneAnnotation.rds")),
    genomeAnnotation = readRDS(paste0(anno_path,"genomeAnnotation.rds")),
    copyArrows = FALSE
  )
  getAvailableMatrices(proj)
  
  ## get gene score matrix
  tmp=getMatrixFromProject(proj,useMatrix = "GeneScoreMatrix")
  genematrix= tmp@assays@data[[1]] #dgcMatrix
  rownames(genematrix)= tmp@elementMetadata$name #genename of dgcMatrix
  rm(tmp); dim(genematrix)
  write.table(genematrix,paste0("~/scATAC/reslt_scATAC/clustering/",sample.name[i],"_ArchR_geneActivity.txt"), col.names = TRUE, row.names = TRUE, sep = "\t")
  
  ## doublet assignment
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
  )
  projHeme2 <- filterDoublets(proj)
  doublets = rownames(getCellColData(proj))
  doublets = doublets[!(doublets %in% rownames(getCellColData(projHeme2)))]
  write.table(doublets,paste0("~/scATAC/reslt_scATAC/clustering/",sample.name[i],"_ArchR_doublets.txt"), col.names = TRUE, row.names = TRUE, sep = "\t")
  
  ## save
  saveArchRProject(ArchRProj = proj,outputDirectory=paste0("~/scATAC/reslt_scATAC/clustering/",sample.name[i]))
}

sessionInfo()

