library(Signac); library(Seurat)

samples=list()
samples[["limb"]]=c(1,2);samples[["nasal"]]=c(3,4);samples[["somite"]]=c(5,6)
sample.name <- c("L21","L24","N15","N18","S12","S15")
data.name <- c("10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_10042020","10x_atac_12012021","10x_atac_12012021")
anno_path <- "/scicore/home/tschoppp/GROUP/references/genomes/ENS_g6/chicken_Gg6_atac_pre/"
rds_path <- paste0("~/scATAC/integration/",paste0(names(samples),"_atac_integrated010721_meta.rds"))

for(i in 1:length(samples)){
  ## read in meta file
  tmp.samples=names(samples)[i]
  hm.integrated=readRDS(rds_path[i])
  hm.integrated$dataset=as.character(hm.integrated$dataset)
  hm.integrated$fine=as.character(hm.integrated$fine)
  for(m in unique(hm.integrated$dataset)){
    for(n in unique(hm.integrated$fine)){
      stage=m
      cls=n
      cells = rownames(hm.integrated[hm.integrated$dataset == stage & hm.integrated$fine == cls,])
      cells = gsub(paste0(stage,"_"),"",cells)
      write.table(cells,paste0("~/scATAC/peakcall/",tmp.samples,"/",stage,"_cls",cls,".txt"),col.names = F,row.names = F,quote = F)
    }
  }
  
  ## subset-bam
  # https://github.com/10XGenomics/subset-bam
  data_path <- paste0("/scicore/home/tschoppp/GROUP/mapped_data/",data.name[i],"/",sample.name[i],"/outs/")
  for(m in unique(hm.integrated$dataset)){
    for(n in unique(hm.integrated$fine)){
      stage=m
      cls=n
      input=paste0(data_path,"possorted_bam.bam")
      cid=paste0("~/scATAC/peakcall/",tmp.samples,"/",stage,"_cls",cls,".txt");print(paste0("run with: ",cid))
      output=paste0("~/scATAC/peakcall/",tmp.samples,"/",stage,"_cls",cls,".bam")
      cmd=paste0("sbatch ~/Rscript/testGit/scATAC/peakcall_subsetBam.sh -b ",input," -c ",cid," -o ",output)
      system(cmd)
    }
  }
  
  ## merge bam of 2 stages
  #samtools merge merged.bam bam1.bam bam2.bam
  for(n in unique(hm.integrated$fine)){
    cls=n
    index=c("-m","-n")
    x=1;input=""
    for(m in unique(hm.integrated$dataset)){
      tmp=paste0("~/scATAC/peakcall/",tmp.samples,"/",m,"_cls",cls,".bam")
      input=paste0(input," ",index[x]," ",tmp);x=x+1
    }
    output=paste0("~/scATAC/peakcall/",tmp.samples,"/",tmp.samples,"_cls",cls,".bam")
    cmd=paste0("sbatch ~/Rscript/testGit/scATAC/peakcall_mergebam.sh -o ",output,input)
    system(cmd)
  }
  
  ## macs2 callpeak
  # samtools index $bamfile
  # macs2 callpeak -t $bamfile -f BAMPE -n $output -g 1.06e9 --keep-dup all --nomodel --shift 100 --extsize 200 --call-summits
  for(n in unique(hm.integrated$fine)){
    cls=n
    input=paste0("~/scATAC/peakcall/",tmp.samples,"/",tmp.samples,"_cls",cls,".bam");print(paste0("run with: ",input))
    output=paste0("~/scATAC/peakcall/",tmp.samples,"/",tmp.samples,"_cls",cls,".peak")
    cmd=paste0("sbatch ~/Rscript/testGit/scATAC/peakcall_macs2.sh -b ",input," -o ",output)
    system(cmd)
  }
  
}

sessionInfo()

