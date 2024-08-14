file_paths <- list.files("/media/baranov_lab/New Volume/Cornea_Seurat_Objs/", pattern = "*.h5Seurat", full.names = TRUE)
seurat_objects <- list()
file_paths[25]
data1 <- LoadH5Seurat(file_paths[1])
data2 <- LoadH5Seurat(file_paths[2])
data3 <- LoadH5Seurat(file_paths[3])
data4 <- LoadH5Seurat(file_paths[4])
data5 <- LoadH5Seurat(file_paths[5])
data6 <- LoadH5Seurat(file_paths[6])
data7 <- LoadH5Seurat(file_paths[7])
data8 <- LoadH5Seurat(file_paths[8])
data9 <- LoadH5Seurat(file_paths[9])
data10 <- LoadH5Seurat(file_paths[10])
data11 <- LoadH5Seurat(file_paths[11])
data12 <- LoadH5Seurat(file_paths[12])
data13 <- LoadH5Seurat(file_paths[13])
data14 <- LoadH5Seurat(file_paths[14])
data15 <- LoadH5Seurat(file_paths[15])
data16 <- LoadH5Seurat(file_paths[16])
data17 <- LoadH5Seurat(file_paths[17])
data18 <- LoadH5Seurat(file_paths[18])
data19 <- LoadH5Seurat(file_paths[19])
data20 <- LoadH5Seurat(file_paths[20])
data21 <- LoadH5Seurat(file_paths[21])
data22 <- LoadH5Seurat(file_paths[22])
data23 <- LoadH5Seurat(file_paths[23])
data24 <- LoadH5Seurat(file_paths[24])
data25 <- LoadH5Seurat(file_paths[25])
data26 <- LoadH5Seurat(file_paths[26])
data27 <- LoadH5Seurat(file_paths[27])
data28 <- LoadH5Seurat(file_paths[28])
data29 <- LoadH5Seurat(file_paths[29])
data30 <- LoadH5Seurat(file_paths[30])
data31 <- LoadH5Seurat(file_paths[31])
data32 <- LoadH5Seurat(file_paths[32])
data33 <- LoadH5Seurat(file_paths[33])
data34 <- LoadH5Seurat(file_paths[34])
data35 <- LoadH5Seurat(file_paths[35])
data36 <- LoadH5Seurat(file_paths[36])
data37 <- LoadH5Seurat(file_paths[37])
data38 <- LoadH5Seurat(file_paths[38])
data39 <- LoadH5Seurat(file_paths[39])
data40 <- LoadH5Seurat(file_paths[40])
data41 <- LoadH5Seurat(file_paths[41])
gc()
data1$batch <- 'Batch41'
data2$batch <- 'Batch18'
data3$batch <- 'Batch19'
data4$batch <- 'Batch20'
data5$batch <- 'Batch29'
data6$batch <- 'Batch30'
data7$batch <- 'Batch21'
data8$batch <- 'Batch22'
data9$batch <- 'Batch23'
data10$batch <- 'Batch24'
data11$batch <- 'Batch25'
data12$batch <- 'Batch28'
data13$batch <- 'Batch26'
data14$batch <- 'Batch27'
data15$batch <- 'Batch31'
data16$batch <- 'Batch32'
data17$batch <- 'Batch33'
data18$batch <- 'Batch34'
data19$batch <- 'Batch35'
data20$batch <- 'Batch36'
data21$batch <- 'Batch17'
data22$batch <- 'Batch16'
data23$batch <- 'Batch15'
data24$batch <- 'Batch14'
data25$batch <- 'Batch13'
data26$batch <- 'Batch12'
data27$batch <- 'Batch37'
data28$batch <- 'Batch38'
data29$batch <- 'Batch39'
data30$batch <- 'Batch1'
data31$batch <- 'Batch2'
data32$batch <- 'Batch3'
data33$batch <- 'Batch4'
data34$batch <- 'Batch5'
data35$batch <- 'Batch6'
data36$batch <- 'Batch7'
data37$batch <- 'Batch8'
data38$batch <- 'Batch9'
data39$batch <- 'Batch10'
data40$batch <- 'Batch11'
data41$batch <- 'Batch40'
integration_list <- list(
  data1, data2, data3, data4, data5, data6, data7, data8, data9, data10,
  data11, data12, data13, data14, data15, data16, data17, data18, data19, data20,
  data21, data22, data23, data24, data25, data26, data27, data28, data29, data30,
  data31, data32, data33, data34, data35, data36, data37, data38, data39, data40,
  data41)
# Verify the integration list
print(length(integration_list))  # Should print 102
print(names(integration_list))   # Should print NULL since we didn't name the elements
gc()
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 2000)
gc()
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.anchors
data.anchorsrpca <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features, reduction = 'rpca')
data.anchorsrpca
gc()
cornea <- IntegrateData(anchorset = data.anchors)
gc()
cornea <- ProcessInt(cornea)
ProcessInt <- function(data.integrated){
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}
DefaultAssay(cornea) <- 'integrated'
options(future.globals.maxSize = 100000 * 1024^2)
cornea150 <- ProcessInt(cornea)
markerscornea <- FindAllMarkers(cornea, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
gc()
markerscornea150 <- FindAllMarkers(cornea150, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
DimPlot(cornea150, label = T, repel = T, label.box = T, raster = T)
test <- markerscornea150 %>%
  group_by(cluster) %>%
  slice_max(n=25, order_by = avg_log2FC)
write.csv(test, '/home/baranov_lab/10X/cornea150_markers.csv')
DefaultAssay(cornea) <- 'RNA'
DefaultAssay(cornea150) <- 'RNA'
DotPlot(cornea150, features = toupper(c('Rho',
                             'Rcvrn',
                             'Pde6b',
                             'Arr3',
                             'Pde6h',
                             'Opn1sw',
                             'Grm6',
                             'Grik1',
                             'Grik2',
                             'Pcp2',
                             'Pcp4',
                             'Prkca',
                             'Vsx2',
                             'Vsx1',
                             'Cck',
                             'Otx2',
                             'Gad1',
                             'Gad2',
                             'Rbpms',
                             'Pou4f1',
                             'Pou4f2',
                             'Pou4f3',
                             'Sncg',
                             'Nefl',
                             'Thy1',
                             'Tfap2a',
                             'Tfap2b',
                             'Slc32a1',
                             'Onecut1',
                             'Onecut2',
                             'Sln',
                             'Rlbp1',
                             'Rpe65',
                             'Gfap',
                             'Sox2',
                             'Pax6',
                             'Slc1a3',
                             'Cx3cr1',
                             'Cd74',
                             'Emcn',
                             'Pdgfra',
                             'Pdgfrb',
                             'Vcan',
                             "Sox10",
                             'Mog',
                             'Myrf',
                             'Mag',
                             'Plp1',
                             'Olig2')))
DotPlot(cornea150, features = toupper(c('Sox8',
                             'Glt1',
                             'Glast',
                             'Mbp',
                             'Olig1',
                             'Ng2',
                             'Glul',
                             'Aqp4','Cd8b', 'Cdba', 'Il32', 'Cd4', 'Il7r', 'Cd3g', 'Cd3d',
                             'Basp', 'Cxcr2', 'Nampt', 'Fcgr3b', 'Csf3r', 'Ms4a2', 'Kit', 'Gata2', 'Cpa3', 'Ilirl1', 'Tpsd1', 'Tpsab1', 'Tpsb2', 'Ahr', 'Klrb1', 'Mtrnr2l1', 'Mtrnr2l10', 'Mtrnr2l6', 'Fn1', 'C1qb', 'Apoc', 'Apoe', 'S100a3', 'Aif1', 'Cd14', 'Lyz', 'Fcn1', 'Vcan', 'Gzmb', 'Gzma', 'Gzmh', 'Cst7', 'Klrf1', 'Fgfbp2', 'Nkg7', 'Gnly', 'Ilsra', 'Fceria', 'Csf2ra', 'Cpvl', 'Gsts', 'Cd1c', 'Amfr', 'Cd22', 'Ralgps2', 'Tcl1abank1', 'Ms4a1', 'Cd79b', 'Cd79a')))
library(SeuratDisk)
SaveH5Seurat(cornea150, '/media/baranov_lab/New Volume/cornea150.h5Seurat', overwrite = TRUE)
SaveH5Seurat(cornea, '/media/baranov_lab/New Volume/cornea30.h5Seurat', overwrite = TRUE)
saveRDS(cornea150, '/media/baranov_lab/New Volume/cornea150.rds')
saveRDS(cornea, '/media/baranov_lab/New Volume/cornea30.rds')
