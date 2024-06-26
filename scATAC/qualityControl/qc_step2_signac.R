library(Signac); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)

sample.name <- c("L21","L24","N15","N18","S12","S15")
data.name <- c("10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_12012021")
anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/chicken_Gg6_atac_pre/"

library(BSgenome.Ggallus.ensembl.galGal6)
genome=seqlengths(BSgenome.Ggallus.ensembl.galGal6)

db <- ensDbFromGtf(gtf=paste0(anno_path,"Gg6_extended_200819.filter.gtf"), organism = "Gallus_gallus", genomeVersion = "GRCg6a", version = 97)
edb <- EnsDb(db); rm(db)
gene.ranges <- genes(edb); length(gene.ranges$gene_id)
annotations <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]

#remove MT, W chromosome and other small contigs
#dropSeqlevels(gr, "W")
genome=genome[names(genome) %in% c(1:33,"Z")]

for(i in 1:length(sample.name)){
  data_path <- paste0("/scicore/home/tschoppp/GROUP/mapped_data/",data.name[i],"/",sample.name[i],"/outs/")
  
  ## input to Signac
  counts <- Read10X_h5(filename = paste0(data_path,"filtered_peak_bc_matrix.h5")); dim(counts)
  metadata <- read.csv(file = paste0(data_path,"singlecell.csv"),header = TRUE,row.names = 1)
  chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = paste0(data_path,'fragments.tsv.gz'), min.cells = 10, min.features = 200)
  
  mtx <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  Annotation(mtx) <- annotations
  
  ## Computing QC Metrics
  mtx <- NucleosomeSignal(object = mtx)
  mtx <- TSSEnrichment(object = mtx, fast = FALSE)
  mtx$pct_reads_in_peaks <- mtx$peak_region_fragments / mtx$passed_filters * 100
  
  mtx$high.tss <- ifelse(mtx$TSS.enrichment > 2, 'High', 'Low');table(mtx$high.tss)
  p1=TSSPlot(mtx, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
  plot(p1)
  
  mtx$nucleosome_group <- ifelse(mtx$nucleosome_signal > 4, 'NS > 4', 'NS < 4');table(mtx$nucleosome_group)
  p2=FragmentHistogram(object = mtx, group.by = 'nucleosome_group',region = "1-1-190000000")
  plot(p2)
  
  p3=VlnPlot(object = mtx,
          features = c('pct_reads_in_peaks', 'peak_region_fragments',
                       'nucleosome_signal'), #, 'TSS.enrichment'
          pt.size = 0.1,
          ncol = 3) + NoLegend()
  plot(p3)
  
  ## Finally we remove cells that are outliers for these QC metrics.
  mtx.fl <- subset(mtx, subset = peak_region_fragments > 1000 & 
                     peak_region_fragments < 100000 & 
                     pct_reads_in_peaks > 15 &
                     nucleosome_signal < 4 & TSS.enrichment > 2)
  
  ## add information of doublets (from ArchR)
  doublets=read.table(paste0("/scicore/home/tschoppp/wang0007/scATAC/reslt_scATAC/clustering/",sample.name[i],"_ArchR_doublets.txt"),
                      header = TRUE, row.names = 1, sep = "\t")
  doublets=as.character(doublets[,1])
  doublets=gsub(paste0(sample.name[i],"#"),"",doublets)
  mtx.fl@meta.data$doublet = "F"
  mtx.fl@meta.data[rownames(mtx.fl@meta.data) %in% doublets,"doublet"]="T"
  table(mtx.fl$doublet)
  mtx.fl = subset(mtx.fl, subset = doublet == "F"); mtx.fl
  
  print(paste0("the number of cells before QC:",ncol(mtx)))
  print(paste0("the number of cells after QC:",ncol(mtx.fl)))
  
  ## add 5-kb bin matrix
  binmatrix=GenomeBinMatrix(
    fragments = Fragments(mtx.fl),
    genome = genome,
    cells = colnames(mtx.fl),
    binsize = 5000, verbose = F
  )
  mtx.fl[['bins']] <- CreateAssayObject(counts = binmatrix)
  rm(binmatrix)
  
  ## save
  saveRDS(mtx.fl,paste0("~/scATAC/reslt_scATAC/clustering/",sample.name[i],"_atac_",format(Sys.Date(), "%d%m%y"),".rds"))
  
}

sessionInfo()

