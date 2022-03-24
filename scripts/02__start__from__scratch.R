
# load seurat object ------------------------------------------------------

load("/workdir/sct__harmony__split__by__orig.ident.Rdata")

sobj.sct <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)

VariableFeatures(sobj.sct) <- var.features

sobj.sct <- RunPCA(sobj.sct, verbose = FALSE, dims = 1:50)
sobj.sct = RunHarmony.Seurat(sobj.sct, 'orig.ident', assay.use = "SCT", 
                             project.dim = FALSE, verbose = TRUE, max.iter.harmony = 20, plot_convergence = T)
sobj.sct = FindNeighbors(sobj.sct, reduction = "harmony", dims = 1:50)
sobj.sct = FindClusters(sobj.sct, resolution = 0.6, graph.name = "SCT_snn")
sobj.sct[['umap']] <- RunUMAP2(Embeddings(sobj.sct, 'harmony')[, 1:50], 
                               assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)

Idents(sobj.sct) <- sobj.sct$SCT_snn_res.0.6
DimPlot(sobj.sct, label = T, reduction = "umap")

saveRDS(sobj.sct, "data/sobj.sct.Rds")

# create symphony ref object ----------------------------------------------

ref <- buildReferenceFromSeurat(sobj.sct, assay = 'SCT', 
                                verbose=TRUE, save_umap=TRUE, save_uwot_path='cache_symphony_sct.uwot')
saveRDS(ref, "data/ref-sv4-res0.6.RDS")


query <- mapQuery(
  q@assays$SCT@scale.data, 
  q@meta.data, 
  ref.sv4.res0.6,
  vars = 'orig.ident', 
  do_normalize = FALSE,
  return_type = 'Seurat' # return a Seurat object
)

options(repr.plot.height = 4, repr.plot.width = 10)
(DimPlot(ref.sv4.res0.6, reduction = 'umap', shuffle = TRUE) + labs(title = 'Original Reference (Clusters)')) + 
  (DimPlot(q, reduction = 'umap', shuffle = TRUE) + labs(title = 'Mapped Query (Donors)'))



q <- knnPredict.Seurat(q, ref.sv4.res0.6, 'SCT_snn_res.0.6')
DimPlot(q, reduction = 'umap', group.by = 'SCT_snn_res.0.6', shuffle = TRUE)
