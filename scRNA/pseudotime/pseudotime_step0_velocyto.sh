# do the same for all 9 samples
velocyto run -@ 8 --samtools-memory 4000 \
  -e L21 -o /scicore/home/tschoppp/wang0007/scRNA/loom \
  -b /scicore/home/tschoppp/wang0007/scRNA/loom/L21_barcodes.tsv \
  /scicore/home/tschoppp/GROUP/mapped_data/10x_L21_ENSG6_extended_260819/outs/possorted_genome_bam.bam \
  /scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/Gg6_extended_200819.gtf
