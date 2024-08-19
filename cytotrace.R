cytotrace2_result_fullmodel <- cytotrace2(Ziming1,
                                          species = "human",
                                          is_seurat = TRUE,
                                          slot_type = "counts",
                                          full_model = TRUE,
                                          batch_size = 10000,
                                          smooth_batch_size = 1000,
                                          parallelize_models = TRUE,
                                          parallelize_smoothing = TRUE,
                                          ncores = 32,
                                          max_pcs = 200,
                                          seed = 14)

annotation <- data.frame(phenotype = Ziming1@meta.data$EK_PB_anno_2024) %>% set_rownames(., colnames(Ziming1))
# plotting
plots5 <- plotData(cytotrace2_result = cytotrace2_result_fullmodel,
                   annotation = annotation,
                   is_seurat = TRUE)
plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- Ziming1@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- Ziming1@reductions[["umap"]]@cell.embeddings[,2]

plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- Ziming1@reductions[["umap"]]@cell.embeddings[,1]
plots5[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- Ziming1@reductions[["umap"]]@cell.embeddings[,2]



plots5$CytoTRACE2_UMAP
plots5$CytoTRACE2_Potency_UMAP
plots5$Phenotype_UMAP
plots5$CytoTRACE2_Boxplot_byPheno
DimPlot(Ziming1)
