## Normal package import
library(dplyr)
library(patchwork)
library(ggcorrplot)
library(reshape2)
library(ggplot2)

## Seurat based
library(Seurat)

## contamination
library(decontX)

## Doublet
library(DoubletFinder)

## GO/KEGG
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Rgraphviz)

## Annotation
library(SingleR)


pipe_qc = function(mydata){
  ### 01_QC ###-------------------------------------------------------------------------------------------------------------------------------
  if (file.exists("01_QC") == F){
    dir.create("01_QC")
  }else { }
  setwd("01_QC")

  mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")
  mydata[["percent.ribo"]] <- PercentageFeatureSet(mydata, pattern = "^RP[SL]")
  mydata[["percent.hb"]] <- PercentageFeatureSet(mydata, pattern = "^HB[^(P)]")

  vlnplot=VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"), ncol = 5)
  ggsave(filename = "qc_vlnplot.pdf", plot = vlnplot,device = 'pdf',dpi = 300, width = 15, height = 8)
  ggsave(filename = "qc_vlnplot.png", plot = vlnplot,device = 'png',dpi = 300, width = 15, height = 8)

  print('01_QC process finished')
  setwd("..")

  return(mydata)
}


pipe_clustering = function(mydata){
  ### 02_cluster ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("02_clustering") == F){
    dir.create("02_clustering")
  }else { }
  setwd("02_clustering")

  mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000, verbose = F)
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- FindNeighbors(mydata, dims = 1:15, verbose = F)
  mydata <- FindClusters(mydata, resolution = 0.1, verbose = F)
  mydata <- RunUMAP(mydata, dims = 1:15, verbose = F)
  umap=DimPlot(mydata, reduction = "umap",group.by = "seurat_clusters",label=T)
  ggsave(filename = "clusters_umap.pdf", plot = umap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_umap.png", plot = umap,device = 'png',dpi = 300, width = 9, height = 8)

  cell_clustering=cbind(colnames(mydata),mydata@meta.data$seurat_clusters)
  write.csv(cell_clustering,file = "cell_clustering.csv")

  # pheatmap
  av<-AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cor_heatmap = pheatmap::pheatmap(cor(av[cg,]))
  ggsave(filename = "clusters_cor_heatmap.pdf", plot = cor_heatmap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_cor_heatmap.png", plot = cor_heatmap,device = 'png',dpi = 300, width = 9, height = 8)

  print('02_clustering process finished')
  setwd("..")

  return(mydata)
}


pipe_cellstatus = function(mydata){
  ### 03_check abnormal cell ###------------------------------------------------------------------------------------------------------------
  if (file.exists("03_cellstatus") == F){
    dir.create("03_cellstatus")
  }else { }
  setwd("03_cellstatus")

  ###### 03_sub contaminated RNA ######
  counts <- mydata@assays$RNA@counts
  decontX_results <- decontX(counts)
  mydata$Contamination = decontX_results$contamination
  mydata$Contam = ifelse(decontX_results$contamination<0.25,'not_contaminated','contaminated')
  Contamplot <- FeaturePlot(mydata,features = 'Contamination', pt.size = .1, reduction = 'umap')
  ggsave(filename = "contamination_featureplot.pdf", plot = Contamplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "contamination_featureplot.png", plot = Contamplot,device = 'png',dpi = 300, width = 9, height = 8)
  Contamplot <- DimPlot(mydata, pt.size = .1, reduction = 'umap', group.by = "Contam")
  ggsave(filename = "contamination_umap.pdf", plot = Contamplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "contamination_umap.png", plot = Contamplot,device = 'png',dpi = 300, width = 9, height = 8)

  ###### 03_sub DoubletFinder ######
  mpK <- 10
  annotations <- mydata@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.04*ncol(mydata@assays$RNA@data))
  seurat_filterDouble <- doubletFinder_v3(mydata, PCs = 1:15, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  cn<-colnames(seurat_filterDouble@meta.data)
  cn[length(cn)] <- "Doublet"
  colnames(seurat_filterDouble@meta.data)<-cn
  mydata$Doublet=seurat_filterDouble$Doublet
  Doubletplot <- DimPlot(mydata, pt.size = .1, reduction = 'umap', group.by = "Doublet")
  ggsave(filename = "doublet_umap.pdf", plot = Doubletplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "doublet_umap.png", plot = Doubletplot,device = 'png',dpi = 300, width = 9, height = 8)

  ###### 03_sub cell cycling ######
  mydata<- CellCycleScoring(mydata, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  cyclingplot <- DimPlot(mydata, pt.size = .1, reduction = 'umap', group.by = "Phase")
  ggsave(filename = "cellcycling_umap.pdf", plot = cyclingplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "cellcycling_umap.png", plot = cyclingplot,device = 'png',dpi = 300, width = 9, height = 8)

  print('03_cellstatus process finished')
  setwd("..")

  return(mydata)
}


pipe_DEGs = function(mydata, by = 'seurat_clusters'){
  ### 04_DEG ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("04_DEGs") == F){
    dir.create("04_DEGs")
  }else { }
  setwd("04_DEGs")

  ###### 04_sub findmarkers ######
  mydata@active.ident <- as.factor(mydata$seurat_clusters)
  pbmc.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.01, verbose = FALSE)
  Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  degs_heatmap<-DoHeatmap(mydata, features = Top5$gene, angle = -50, hjust=0.8)
  ggsave(filename = "TopDEGs_heatmap.pdf", plot = degs_heatmap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "TopDEGs_heatmap.png", plot = degs_heatmap,device = 'png',dpi = 300, width = 9, height = 8)
  write.csv(pbmc.markers, file='clusters_DEGs.csv')

  ###### 04_sub featureplot ######
  Top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  for (i in levels(Top20$cluster)){
    name <- paste("cluster",i,sep="")
    dir.create(name)
    sub <- subset(Top20, cluster %in% i)
    for(j in sub$gene){
      p<-FeaturePlot(mydata,features=j,pt.size = .1, reduction = 'umap')
      ggsave(filename = paste(name,"/",j,".pdf",sep=""), plot = p,device = 'pdf',dpi = 300, width = 9, height = 8)
      ggsave(filename = paste(name,"/",j,".png",sep=""), plot = p,device = 'png',dpi = 300, width = 9, height = 8)
    }
  }

  print('04_DEGs process finished')
  setwd("..")

  return(mydata)
}


pipe_pathway = function(mydata, by = 'seurat_clusters'){
  ### 05_GOKEGG ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("05_GOKEGG") == F){
    dir.create("05_GOKEGG")
  }else { }
  setwd("05_GOKEGG")

  if (file.exists("../04_DEGs/clusters_DEGs.csv") == T){
    group_g <- read.csv('../04_DEGs/clusters_DEGs.csv', header = T, row.names = 1)
  }else {
    mydata@active.ident <- as.factor(mydata$seurat_clusters)
    group_g <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.01, verbose = FALSE)
  }
  tmp <- bitr(group_g$gene, fromType="SYMBOL",
              toType = c("ENSEMBL", "ENTREZID"),
              OrgDb="org.Hs.eg.db")
  de_gene_clusters=merge(tmp,group_g,by.x='SYMBOL',by.y='gene')

  ###### 05_sub GO enrichment ######
  formula_res <- compareCluster(
    ENTREZID~cluster,
    data=de_gene_clusters,
    fun="enrichGO",
    OrgDb="org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)

  t<-formula_res@compareClusterResult
  write.csv(t,file="formula_res_go.csv")

  lineage1_ego <- clusterProfiler::simplify(
    formula_res,
    cutoff=0.1,
    by="p.adjust",
    select_fun=min
  )

  formula_res %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = p.adjust) -> formula_res2

  p = enrichplot::dotplot(formula_res2, showCategory=3)
  ggsave(filename = 'enrichplot_go.pdf', plot = p,device = 'pdf',dpi = 300, width = 12, height = 11)
  ggsave(filename = 'enrichplot_go.png', plot = p,device = 'png',dpi = 300, width = 12, height = 11)


  ###### 05_sub KEGG enrichment ######
  options(clusterProfiler.download.method = "wget")
  formula_res <- compareCluster(
    ENTREZID~cluster,
    data=de_gene_clusters,
    fun="enrichKEGG",
    organism="hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  formula_res %>%
    group_by(cluster) %>%
    top_n(n = 3, wt = p.adjust) -> formula_res2

  t<-formula_res@compareClusterResult
  write.csv(t,file="formula_res_kegg.csv")

  p = enrichplot::dotplot(formula_res2, showCategory=3)
  ggsave(filename = 'enrichplot_kegg.pdf', plot = p,device = 'pdf',dpi = 300, width = 12, height = 11)
  ggsave(filename = 'enrichplot_kegg.png', plot = p,device = 'png',dpi = 300, width = 12, height = 11)

  ###### 05_sub forloop enrichment ######
  options(clusterProfiler.download.method = "wget")
  for (i in levels(mydata$seurat_clusters)){
    dir_name=paste("cluster",i,sep="")
    dir.create(dir_name)
    markers <- subset(pbmc.markers, cluster==i)$gene
    genelist = bitr(markers, fromType="SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    ego <- enrichGO(genelist$ENTREZID,OrgDb = "org.Hs.eg.db", ont = "BP",readable = T,pvalueCutoff = 0.5, qvalueCutoff = 1)
    df_ego <- summary(ego)
    if(dim(df_ego)[1]!=0){
      p1=barplot(ego,title="enrichGO")
      ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".pdf",sep=""), plot = p1,device = 'pdf',dpi = 300, width = 12, height = 11)
      ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".png",sep=""), plot = p1,device = 'png',dpi = 300, width = 12, height = 11)
      t1=ego@result
      write.csv(t1,file=paste(dir_name,"/DEG_GO_cluster",i,".csv",sep=""))
    }
    kk <- enrichKEGG(gene = genelist$ENTREZID, organism ="hsa", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1, use_internal_data =FALSE)
    df_kk <- summary(kk)
    if(dim(df_kk)[1]!=0){
      p2=barplot(kk,title="KEGG")
      ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".pdf",sep=""), plot = p2,device = 'pdf',dpi = 300, width = 12, height = 11)
      ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".png",sep=""), plot = p2,device = 'png',dpi = 300, width = 12, height = 11)
      t2=kk@result
      write.csv(t2,file=paste(dir_name,"/DEG_KEGG_cluster",i,".csv",sep=""))
    }
  }
  print('05_GOKEGG process finished')
  setwd("..")

  return(mydata)

}


advanced_featureplot = function(mydata, genelist, reduction='umap', ncol=4){

  GeneExp <- FetchData(mydata,vars = genelist)
  pc <- Embeddings(mydata,reduction = reduction) %>% data.frame()
  colnames(pc) <- c('Dim_1','Dim_2')
  gbid <- cbind(pc,GeneExp)
  gbidlong <- melt(gbid,id.vars = c('Dim_1','Dim_2'),value.name = 'exp',variable.name = 'gene')

  ggplot(gbidlong,aes(x = Dim_1,y = Dim_2,color = exp)) +
    geom_point(size = 0,show.legend = T) +
    scale_color_gradient(low = 'lightgrey',high = '#DD0000',name = 'Expr') +
    theme_bw(base_size = 24) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          strip.background = element_rect(colour = NA,fill = NA),
          axis.text = element_blank(),
          plot.title = element_text(size=24,hjust = 0.5)) +
    facet_wrap(~gene,ncol = ncol)
}


pipe_anno = function(mydata){
  ### 06_Annotation ###----------------------------------------------------------------------------------------------------------------------
  if (file.exists("06_Anno") == F){
    dir.create("06_Anno")
  }else { }
  setwd("06_Anno")

  ###### 06_sub singleR ######
  refdata <- HumanPrimaryCellAtlasData()
  testdata <- GetAssayData(mydata, slot="data")
  clusters <- mydata$seurat_clusters
  cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                      method = "cluster", clusters = clusters,
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

  mydata$SingleR = "NA"
  for(i in 1:nrow(celltype)){
    mydata@meta.data[which(mydata$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
  }

  p<-DimPlot(mydata, group.by="SingleR", label=T, label.size=0.2, reduction='umap')
  ggsave(filename = "SingleR_anno_umap.pdf", plot = p,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "SingleR_anno_umap.png", plot = p,device = 'png',dpi = 300, width = 9, height = 8)

  ###### 06_sub cellMarkers dotplot ######
  Epi=c('EPCAM', 'KRT5','KRT18')
  Myeloid=c('CD14','CD68','FCGR3A','LYZ')
  DC=c('CD1C','FCER1A','CST3')
  Fib=c('COL1A1','COL1A2','DCN')
  Endo=c('PECAM1','CLDN5','RAMP2')
  T=c('CD3D','CD3E','CD4')
  NK=c('NKG7','KLRD1','GZMB')
  B=c('CD19','CD79A','CD79B')
  Mast=c('TPSB2','CPA3','TPSAB1')
  Cycling=c('MKI67','PCNA','TOP2A')
  gene_list = list(
    Epi=Epi,
    Myeloid=Myeloid,
    DC=DC,
    Fib=Fib,
    Endo=Endo,
    T=T,
    NK=NK,
    B=B,
    Mast=Mast,
    Cycling=Cycling
  )

  p = DotPlot(mydata, assay = "RNA", features = gene_list, group.by = 'seurat_clusters') +
    theme(axis.text.x = element_text(angle = 45,  vjust = 0.9, hjust=0.9)) +
    scale_colour_gradient2(low = "steelblue", mid = "lightgrey", high = "#DD0000") +
    RotatedAxis()
  ggsave(filename = "cellMarkers_dotplot.pdf", plot = p,device = 'pdf',dpi = 300, width = 11, height = 9)
  ggsave(filename = "cellMarkers_dotplot.png", plot = p,device = 'png',dpi = 300, width = 11, height = 9)


  ###### 06_sub cellMarkers featureplot ######

  genelist = c('EPCAM','CD14','CD68','COL1A1','PECAM1','CD3D','CD4','CD8','FOXP3','NKG7','CD79A','MS4A1','CPA3','MZBI','MKI67')
  p = advanced_featureplot(mydata, genelist) +
    RotatedAxis()
  ggsave(filename = "cellMarkers_featureplot.pdf", plot = p,device = 'pdf',dpi = 300, width = 11, height = 9)
  ggsave(filename = "cellMarkers_featureplot.png", plot = p,device = 'png',dpi = 300, width = 11, height = 9)



  print('06_Anno process finished')
  setwd("..")

  return(mydata)

}


pipe_save = function(mydata, name){
  ### 07_out ###----------------------------------------------------------------------------------------------------------------------
  if (file.exists("07_out") == F){
    dir.create("07_out")
  }else { }
  setwd("07_out")

  saveRDS(mydata,file=paste(name,'.rds', sep = ''))
  meta<-mydata@meta.data
  write.csv(meta,file=paste(name,'_metadata.csv', sep = ''))

  print('07_out process finished')
  setwd("..")

  return(mydata)
}


pipe_samplerun = function(indir,outdir,samplename=NULL,save=TRUE, midsave=FALSE){

  ### 10X matrix read ###-------------------------------------------------------------------------------------------------------------------------------
  setwd(outdir)
  mydata = Read10X(data.dir = indir)
  if(!is.null(samplename)){
    dir.create(samplename)
    setwd(samplename)
  }else{
    samplename =strsplit(indir, '/')[[1]][length(strsplit(indir, '/')[[1]])]
    dir.create(samplename)
    setwd(samplename)
  }
  mydata <- CreateSeuratObject(counts = mydata, projec=samplename, min.cells = 3, min.features = 200)
  mydata = pipe_qc(mydata)
  mydata = pipe_clustering(mydata)
  mydata = pipe_cellstatus(mydata)
  mydata = pipe_DEGs(mydata)
  mydata = pipe_pathway(mydata)
  mydata = pipe_anno(mydata)
  pipe_save(mydata,name = samplename)

  return(mydata)
}
