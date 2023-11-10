## basic packages:
library(dplyr)
library(patchwork)
library(reshape2)
## visualizing packages:
library(ggcorrplot)
library(ggplot2)
library(pheatmap)
library(cols4all)
## Seurat based
library(Seurat)
## contamination predicting
library(decontX)
## Doublet predicting
library(DoubletFinder)
## pathway analysis
options(connectionObserver = NULL)
library(AnnotationHub)
library(org.Hs.eg.db)
library(clusterProfiler)
library(Rgraphviz)
## celltype predicting
library(SingleR)

pipe_01qc = function(mydata){
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


pipe_02clustering = function(mydata, hvg=2000, pca=15, res=0.1){
  ### 02_cluster ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("02_clustering") == F){
    dir.create("02_clustering")
  }else { }
  setwd("02_clustering")

  mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- FindNeighbors(mydata, dims = 1:pca, verbose = F)
  mydata <- FindClusters(mydata, resolution = res, verbose = F)
  mydata <- RunUMAP(mydata, dims = 1:pca, verbose = F)
  umap=advanced_dimplot(mydata)
  ggsave(filename = "clusters_umap.pdf", plot = umap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_umap.png", plot = umap,device = 'png',dpi = 300, width = 9, height = 8)

  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")

  # pheatmap
  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cor_heatmap = pheatmap(cor(as.matrix(av[cg,])))
  ggsave(filename = "clusters_cor_heatmap.pdf", plot = cor_heatmap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_cor_heatmap.png", plot = cor_heatmap,device = 'png',dpi = 300, width = 9, height = 8)

  print('02_clustering process finished')
  setwd("..")

  return(mydata)
}

pipe_03clustering_qcfilter = function(mydata, mt = 20, contam = 0.25, singlet=TRUE, hvg=2000, pca=15, res=0.1){
  ### 02_cluster ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("03_clustering_qcfilter") == F){
    dir.create("03_clustering_qcfilter")
  }else { }
  setwd("03_clustering_qcfilter")

  if(!is.na(mt)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = percent.mt < mt)
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of mitochondria, ',cellnum,' cells were screened out.', sep = ''))
  }else{}

  if(!is.na(contam)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = Contamination < contam)
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of contamination, ',cellnum,' cells were screened out.', sep = ''))
  }else{}

  if(isTRUE(singlet)){
    cellnum1 = ncol(mydata)
    mydata <- subset(mydata, subset = Doublet == 'Singlet')
    cellnum2 = ncol(mydata)
    cellnum = cellnum1-cellnum2
    print(paste('During the QC process of Doublet, ',cellnum,' cells were screened out.', sep = ''))
  }else{}

  mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = hvg, verbose = F)
  mydata <- ScaleData(mydata, features = rownames(mydata), verbose = F)
  mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata), verbose = F)
  mydata <- FindNeighbors(mydata, dims = 1:pca, verbose = F)
  mydata <- FindClusters(mydata, resolution = res, verbose = F)
  mydata <- RunUMAP(mydata, dims = 1:pca, verbose = F)
  umap=advanced_dimplot(mydata)
  ggsave(filename = "clusters_umap.pdf", plot = umap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_umap.png", plot = umap,device = 'png',dpi = 300, width = 9, height = 8)

  cell_clustering=as.data.frame(cbind(colnames(mydata),mydata@meta.data$seurat_clusters))
  colnames(cell_clustering) = c('barcodes','seurat_clusters')
  write.csv(cell_clustering,file = "cell_clustering.csv")

  # pheatmap
  av=AverageExpression(mydata, group.by = "seurat_clusters", assays = "RNA", verbose = F)
  av=av[[1]]
  cg=names(tail(sort(apply(av,1,sd)),1000))
  cor_heatmap = pheatmap(cor(as.matrix(av[cg,])))
  ggsave(filename = "clusters_cor_heatmap.pdf", plot = cor_heatmap,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "clusters_cor_heatmap.png", plot = cor_heatmap,device = 'png',dpi = 300, width = 9, height = 8)

  print('03_clustering_qcfilter process finished')
  setwd("..")

  return(mydata)
}


pipe_03cellstatus = function(mydata){
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
  Contamplot <- advanced_dimplot(mydata, group.by='Contam')
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
  Doubletplot <- advanced_dimplot(mydata, group.by='Doublet')
  ggsave(filename = "doublet_umap.pdf", plot = Doubletplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "doublet_umap.png", plot = Doubletplot,device = 'png',dpi = 300, width = 9, height = 8)

  ###### 03_sub cell cycling ######
  mydata<- CellCycleScoring(mydata, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  cyclingplot <- advanced_dimplot(mydata, group.by='Phase')
  ggsave(filename = "cellcycling_umap.pdf", plot = cyclingplot,device = 'pdf',dpi = 300, width = 9, height = 8)
  ggsave(filename = "cellcycling_umap.png", plot = cyclingplot,device = 'png',dpi = 300, width = 9, height = 8)

  print('03_cellstatus process finished')
  setwd("..")

  return(mydata)
}


pipe_04DEGs = function(mydata, by = 'seurat_clusters'){
  ### 04_DEG ###--------------------------------------------------------------------------------------------------------------------------
  if (file.exists("04_DEGs") == F){
    dir.create("04_DEGs")
  }else { }
  setwd("04_DEGs")

  ###### 04_sub findmarkers ######
  mydata@active.ident <- as.factor(mydata$seurat_clusters)
  pbmc.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, return.thresh = 0.01, verbose = FALSE)
  Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  degs_heatmap<-DoHeatmap(mydata, features = Top5$gene, angle = -50, hjust=0.8, raster = FALSE)
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


pipe_05pathway = function(mydata, by = 'seurat_clusters'){
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
  # options(clusterProfiler.download.method = "wget")
  # for (i in levels(mydata$seurat_clusters)){
  #   dir_name=paste("cluster",i,sep="")
  #   dir.create(dir_name)
  #   markers <- subset(group_g, cluster==i)$gene
  #   genelist = bitr(markers, fromType="SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
  #   ego <- enrichGO(genelist$ENTREZID,OrgDb = "org.Hs.eg.db", ont = "BP",readable = T,pvalueCutoff = 0.5, qvalueCutoff = 1)
  #   df_ego <- summary(ego)
  #   if(dim(df_ego)[1]!=0){
  #     p1=barplot(ego,title="enrichGO")
  #     ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".pdf",sep=""), plot = p1,device = 'pdf',dpi = 300, width = 12, height = 11)
  #     ggsave(filename = paste(dir_name,"/DEG_GO_cluster",i,".png",sep=""), plot = p1,device = 'png',dpi = 300, width = 12, height = 11)
  #     t1=ego@result
  #     write.csv(t1,file=paste(dir_name,"/DEG_GO_cluster",i,".csv",sep=""))
  #   }
  #   kk <- enrichKEGG(gene = genelist$ENTREZID, organism ="hsa", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1, use_internal_data =FALSE)
  #   df_kk <- summary(kk)
  #   if(dim(df_kk)[1]!=0){
  #     p2=barplot(kk,title="KEGG")
  #     ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".pdf",sep=""), plot = p2,device = 'pdf',dpi = 300, width = 12, height = 11)
  #     ggsave(filename = paste(dir_name,"/DEG_KEGG_cluster",i,".png",sep=""), plot = p2,device = 'png',dpi = 300, width = 12, height = 11)
  #     t2=kk@result
  #     write.csv(t2,file=paste(dir_name,"/DEG_KEGG_cluster",i,".csv",sep=""))
  #   }
  # }
  print('05_GOKEGG process finished')
  setwd("..")

  return(mydata)

}



pipe_06anno = function(mydata){
  ### 06_Annotation ###----------------------------------------------------------------------------------------------------------------------
  if (file.exists("06_Anno") == F){
    dir.create("06_Anno")
  }else { }
  setwd("06_Anno")

  ###### 06_sub singleR ######
  refdata <- celldex::HumanPrimaryCellAtlasData()
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

  p<-advanced_dimplot(mydata, group.by='SingleR')
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
  p = advanced_featureplot(mydata)
  ggsave(filename = "cellMarkers_featureplot.pdf", plot = p,device = 'pdf',dpi = 300, width = 9, height = 23)
  ggsave(filename = "cellMarkers_featureplot.png", plot = p,device = 'png',dpi = 300, width = 9, height = 23)

  ###### 06_sub cellMarkers violinplot ######
  p = advanced_violinplot(mydata)
  ggsave(filename = "cellMarkers_vlnplot.pdf", plot = p,device = 'pdf',dpi = 300, width = 23, height = 9)
  ggsave(filename = "cellMarkers_vlnplot.png", plot = p,device = 'png',dpi = 300, width = 23, height = 9)



  print('06_Anno process finished')
  setwd("..")

  return(mydata)

}


pipe_07save = function(mydata, name){
  ### 07_out ###----------------------------------------------------------------------------------------------------------------------
  if (file.exists("07_out") == F){
    dir.create("07_out")
  }else { }
  setwd("07_out")

  saveRDS(mydata,file=paste(name,'.rds', sep = ''))
  meta<-mydata@meta.data
  write.csv(meta,file=paste(name,'_metadata.csv', sep = ''))

  print('07_save process finished')
  setwd("..")

  return(mydata)
}


advanced_featureplot = function(mydata, genelist=NULL, reduction='umap', ncol=3){
  if(is.null(genelist)){
    genelist = c('EPCAM', 'KRT5','KRT18','CD14','CD68','FCGR3A','LYZ',
                 'CD1C','FCER1A','CST3','COL1A1','COL1A2','DCN',
                 'PECAM1','CLDN5','RAMP2','CD3D','CD3E','CD4',
                 'NKG7','KLRD1','GZMB','CD19','CD79A','CD79B',
                 'TPSB2','CPA3','TPSAB1','MKI67','PCNA','TOP2A')
    genelist = intersect(genelist, rownames(mydata))
  }else{}

  GeneExp <- FetchData(mydata,vars = genelist)
  pc <- Embeddings(mydata,reduction = reduction) %>% data.frame()
  colnames(pc) <- c('Dim_1','Dim_2')
  gbid <- cbind(pc,GeneExp)
  gbidlong <- melt(gbid,id.vars = c('Dim_1','Dim_2'),value.name = 'exp',variable.name = 'gene')

  ggplot(gbidlong,aes(x = Dim_1,y = Dim_2,color = exp)) +
    geom_point(size = 0,show.legend = T) +
    scale_color_gradient(low = 'lightgrey',high = '#DD0000',name = 'Expr') +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1,
          strip.background = element_rect(colour = NA,fill = NA),
          axis.text = element_blank(),
          plot.title = element_text(size=16,hjust = 0.5)) +
    facet_wrap(~gene,ncol = ncol)
}


advanced_violinplot = function(mydata, genelist=NULL, group.by=NULL,color=NULL){
  if(is.null(genelist)){
    genelist = c('EPCAM', 'KRT5','KRT18','CD14','CD68','FCGR3A','LYZ',
                 'CD1C','FCER1A','CST3','COL1A1','COL1A2','DCN',
                 'PECAM1','CLDN5','RAMP2','CD3D','CD3E','CD4',
                 'NKG7','KLRD1','GZMB','CD19','CD79A','CD79B',
                 'TPSB2','CPA3','TPSAB1','MKI67','PCNA','TOP2A')
    genelist = intersect(genelist, rownames(mydata))
  }else{}

  dd = as.data.frame(t(as.matrix(mydata@assays$RNA@data[genelist,])))
  all(rownames(dd)==rownames(mydata@meta.data))

  if(is.null(group.by)){
    dd$group = mydata$seurat_clusters
  }else{
    dd$group = mydata@meta.data[,group.by]
  }

  dd = melt(dd)
  colnames(dd) = c('cluster','genes', 'expr')

  if(is.null(color)){
    mycol<- c4a('20')
  }else{}

  p = ggplot(data = dd,aes(x = expr, y = cluster, fill = cluster)) +
    geom_violin(scale = 'width',
                draw_quantiles= c(0.25, 0.5, 0.75),
                color= 'black',
                size= 0.45,
                alpha= 0.8) +
    facet_grid(cols = vars(genes), scales = 'free_x')+
    scale_fill_manual(values = mycol) + #填充色修改
    scale_x_continuous(breaks = seq(0, 8, by = 4)) + #x轴刻度显示调整
    theme_bw()+
    theme(
      panel.grid = element_blank(), #移除背景网格线
      axis.text.x = element_blank(), #x轴标签大小调整
      axis.text.y = element_text(size = 16), #y轴标签大小调整
      axis.title.x = element_text(size = 16), #x轴标题大小调整
      axis.title.y = element_blank(), #移除y轴标题
      axis.ticks.x = element_blank(),
      strip.background = element_blank(), #移除分面外围边框
      strip.text.x = element_text(size = 16, angle = 60), #分面文本标签倾斜60°
      legend.title = element_text(size = 16), #图例标题大小调整
      legend.text = element_text(size = 15) #图例标签大小调整
    ) +
    labs(x = 'Log Normalized Expression')

  return(p)
}


advanced_dimplot = function(mydata, group.by=NULL, reduction=NULL,color=NULL){
  if(is.null(group.by)){
    group.by = 'seurat_clusters'
  }else{}

  if(is.null(reduction)){
    reduction = 'umap'
  }else{}

  reduct = unlist(mydata@reductions[reduction][[1]])
  umap = reduct@cell.embeddings %>%
    as.data.frame() %>%
    cbind(cell_type = mydata@meta.data[,group.by])
  colnames(umap) = c('DIM_1','DIM_2','clusters')

  if(is.null(color)){
    mycol<- c4a("poly.wright25")
  }else{}

  ggplot(umap,aes(x= DIM_1 , y = DIM_2 ,color = clusters)) +
    geom_point(size = 0.1 , alpha =0.8 ) +
    scale_color_manual(values = mycol)+
    #scale_color_discrete_c4a_cat("carto.safe")+
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white")) +
    theme(
      legend.title = element_blank(), #去掉legend.title
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=20), #设置legend标签的大小
      legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5))) + #设置legend中 点的大小
    geom_segment(aes(x = min(DIM_1) , y = min(DIM_2) ,
                     xend = min(DIM_1) +3, yend = min(DIM_2) ),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
    geom_segment(aes(x = min(DIM_1)  , y = min(DIM_2)  ,
                     xend = min(DIM_1) , yend = min(DIM_2) + 3),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(umap$DIM_1) +1.5, y = min(umap$DIM_2) -1, label = "DIM_1",
             color="black",size = 6 ) +
    annotate("text", x = min(umap$DIM_1) -1, y = min(umap$DIM_2) + 1.5, label = "DIM_2",
             color="black",size = 6, angle=90) +
    theme(legend.position = "right",
          legend.text = element_text(size=23))
}


pipe_samplerun = function(indir,outdir,samplename=NULL,deep_qc = TRUE,save=TRUE, midsave=FALSE){

  ###### 10X matrix read ######
  setwd(outdir)
  if(!is.null(samplename)){
    dir.create(samplename)
    setwd(samplename)
  }else{
    samplename =strsplit(indir, '/')[[1]][length(strsplit(indir, '/')[[1]])]
    dir.create(samplename)
    setwd(samplename)
  }
  mydata = Read10X(data.dir = indir)
  mydata <- CreateSeuratObject(counts = mydata, projec=samplename, min.cells = 3, min.features = 200)
  mydata = pipe_01qc(mydata)
  mydata = pipe_02clustering(mydata)
  mydata = pipe_03cellstatus(mydata)
  if(deep_qc==TRUE){
    mydata = pipe_03clustering_qcfilter(mydata)
  }else{}
  mydata = pipe_04DEGs(mydata)
  mydata = pipe_05pathway(mydata)
  mydata = pipe_06anno(mydata)
  pipe_07save(mydata,name = samplename)

  return(mydata)
}
