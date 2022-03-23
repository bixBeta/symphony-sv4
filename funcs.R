
# UMAP Dotplots -----------------------------------------------------------

getUmapDotPlots = function(sobj, res, feats){
  
  title=deparse(substitute(sobj))
  Idents(object = sobj) <- sobj@meta.data[,paste0("SCT_snn_res.", res)]
  uplot = DimPlot(sobj , group.by = paste0("SCT_snn_res.", res), label = T, label.size = 6, repel = T, raster = F) + ggtitle(paste0(title, " UMAP -- res ", res))
  
  dplot = DotPlot(object = sobj, assay="RNA", features=feats, col.min=-1.5,
                  group.by = paste0("SCT_snn_res.", res), cols="Spectral", dot.scale=3, cluster.idents = T, scale.by = 'size') + 
    RotatedAxis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste0(title, " DotPlot -- res ", res))
  
  cplot = uplot + dplot
  
  
  png(filename = paste0("sv4__umap__dotPlot__res__",res, ".png"), width = 3000, height = 900, res = 100)
  print(cplot)
  dev.off()
  
}


getUmapDotPlots(sobj = sv4, res = 0.6, feats = lv.feats)



# sankey ------------------------------------------------------------------

v4raw.v.v4t = as.matrix(table(sv4.raw.t.meta$SCT_snn_res.0.6, v4t.obj.added.clusts.meta$SCT_snn_res.0.8))
colnames(v4raw.v.v4t) = paste0("res0.8_T_cluster", colnames(v4raw.v.v4t))
rownames(v4raw.v.v4t) = paste0("res0.6_SV4_Cluster", rownames(v4raw.v.v4t))

getSankey <- function(matrix, floor){
  
  x = as.data.frame(matrix) 
  colnames(x) = c("source", "target", "value")
  
  x = x %>% filter(value > floor)
  
  nodes = as.data.frame(c(as.character(x$source), as.character(x$target)) %>% unique())
  
  colnames(nodes) = "name"
  
  x$IDsource=match(x$source, nodes$name)-1 
  x$IDtarget=match(x$target, nodes$name)-1
  
  y = sankeyNetwork(Links = x, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=12, nodePadding=20, iterations = 0)
  
  return(y)
  
}
getSankey(matrix = v4raw.v.v4t, floor = 10)



# nde ---------------------------------------------------------------------

nde = function(de.list, padj){
  
  z = list()
  
  for (i in 1:length(de.list)) {
    z[[i]] = length(which(pluck(de.list, i, "p_val_adj") < padj))
    names(z)[[i]] <- names(de.list)[[i]]
    
  }
  
  nDE.de.list = do.call("rbind",z)
  colnames(nDE.de.list) = deparse(substitute(de.list))
  
  
  return(nDE.de.list)
}


x.1 = nde(de.list = caseD1.controlD1.list, padj = 0.05)
x.2 = nde(de.list = caseD2.controlD2.list, padj = 0.05)
x.3 = nde(de.list = controlD1.controlD2.list, padj = 0.05)
x.4 = nde(de.list = caseD1.caseD2.list, padj = 0.05)

df = cbind(x.1,x.2,x.3,x.4)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(reshape2)
nde.gg = melt(df)
colnames(nde.gg) = c("cluster","condition", "nDE")

ggplot(nde.gg, aes(x=cluster, y=nDE , size = nDE , color = condition)) + scale_colour_viridis_d("condition") +
  facet_wrap(~condition) + ggtitle("Number of DE genes - SV4 - Tcells res 0.8") +
  geom_point(alpha=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

