source("symohony-utils-seurat.R")
source("funcs.R")
source("packages.R")

# import seurat object ----------------------------------------------------

sv4 = readRDS("/workdir/SV4__Backup/seurat_object_v4.RDS")
sv4.meta = sv4@meta.data
sanity.meta = sv4@meta.data


# create symphony ref -----------------------------------------------------

## rerun harmony to get the MISC slot filled

sv4.harm = RunHarmony.Seurat(sv4, 'orig.ident', assay.use = "SCT", 
                  project.dim = FALSE, verbose = TRUE, max.iter.harmony = 20, plot_convergence = T)


## Currently, Seurat does not let you cache the umap model for future mapping
## Therefore, please use this custom function to learn a saveable UMAP model

sv4.harm[['umap']] <- RunUMAP2(Embeddings(sv4.harm, 'harmony')[, 1:50], 
                          assay='RNA', verbose=FALSE, umap.method='uwot', return.model=TRUE)

sv4.harm = FindNeighbors(sv4.harm, reduction = "harmony", dims = 1:50)
sv4.harm = FindClusters(sv4.harm, resolution = c(0.4, 0.5, 0.6, 0.8), graph.name = "SCT_snn")

sv4.harm.meta = sv4.harm@meta.data

ref <- buildReferenceFromSeurat(sv4.harm, assay = 'SCT', 
                                verbose=TRUE, save_umap=TRUE, save_uwot_path='cache_symphony_sct.uwot')
saveRDS(ref, "data/symphony-ref-sv4-created-from-sv4.harm.RDS")




