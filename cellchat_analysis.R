library(CellChat) #built on v2.1.0
cellchat <- createCellChat(object = Ziming_vivo , group.by = "anno_cc")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) #or cellchat <- smoothData(cellchat, adj = PPI.human) if version > 2.1.0
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

cellchat_z_vivo <- cellchat
saveRDS(cellchat, file = "C://Bioinf/Ziming/cellchat_vivo.rds")
