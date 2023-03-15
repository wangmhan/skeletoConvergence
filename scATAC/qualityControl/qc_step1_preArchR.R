library(ArchR); packageVersion("ArchR")
library(GenomeInfoDb); library(ensembldb)

set.seed(1234)
addArchRThreads(threads = 1) 

##########################
#Creating a Custom ArchRGenome
##########################

#https://www.archrproject.com/bookdown/getting-set-up.html
library(BSgenome.Ggallus.ensembl.galGal6) #library(BSgenome.Ggallus.UCSC.galGal6)
genomeAnnotation = createGenomeAnnotation(genome = "BSgenome.Ggallus.ensembl.galGal6",filterChr = "MT")

anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/chicken_Gg6_atac_pre/"
db <- ensDbFromGtf(gtf=paste0(anno_path,"Gg6_extended_200819.filter.gtf"), 
                   organism = "Gallus_gallus", genomeVersion = "GRCg6a", version = 97)
edb <- EnsDb(db); rm(db)
edb = addFilter(edb, filter = GeneBiotypeFilter('protein_coding',"=="))
edb = addFilter(edb, filter = SeqNameFilter(c('AADN05001525.1','KZ626819.1','KZ626826.1','KZ626830.1','MT',
                                              'KZ626833.1','KZ626834.1','KZ626835.1','KZ626836.1','Z','W','KZ626839.1'),"!=")) #
gene.ranges <- genes(edb); length(gene.ranges$gene_id)
tss.ranges <- GRanges(seqnames = seqnames(gene.ranges),
                      ranges = IRanges(start = start(gene.ranges), width = 2),
                      strand = strand(gene.ranges)
)
exon.ranges = exons(edb)
geneAnnotation = createGeneAnnotation(
  TSS = tss.ranges, 
  exons = exon.ranges, 
  genes = gene.ranges
)
unique(seqnames(geneAnnotation$genes))

saveRDS(genomeAnnotation, paste0(anno_path,"genomeAnnotation.rds"))
saveRDS(geneAnnotation, paste0(anno_path,"geneAnnotation.rds"))


