
library(escape)

gene.sets1 <- getGeneSets(library = "C5", gene.sets = c('GOBP_MICROGLIA_DIFFERENTIATION',
                                                        'GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                                                        'GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY',
                                                        'GOBP_MICROGLIAL_CELL_MIGRATION',
                                                        'GOBP_MICROGLIAL_CELL_PROLIFERATION',
                                                        'GOBP_REGULATION_OF_MICROGLIA_CELL_ACTIVATION',
                                                        'GOBP_PHAGOCYTOSIS',
                                                        'GOBP_IMMUNE_RESPONSE',
                                                        'GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY',
                                                        'GOBP_ACTIVATION_OF_IMMUNE_RESPONSE',
                                                        'GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE',
                                                        'GPBP_INNATE_IMMUNE_RESPONSE',
                                                        'GOBP_INNATE_IMMUNE_RESPONSE_ACTIVATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY',
                                                        'GOBP_INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNALING_PATHWAY',
                                                        'GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE',
                                                        'GOBP_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE ',
                                                        'GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE ',
                                                        'GOBP_T_CELL_CYTOKINE_PRODUCTION ',
                                                        'GO_BP_ASTROCYTE_ACTIVATION',
                                                        'GOBP_MACROPHAGE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                                                        'GOBP_TRYPTOPHAN_CATABOLIC_PROCESS',
                                                        'GOBP_TRYPTOPHAN_CATABOLIC_PROCESS_TO_KYNURENINE',
                                                        'GOBP_TRYPTOPHAN_METABOLIC_PROCESS',
                                                        'GOBP_TRYPTOPHAN_TRANSPORT',
                                                        'GOBP_SEROTONIN_METABOLIC_PROCESS',
                                                        'GOBP_VITAMIN_B6_METABOLIC_PROCESS',
                                                        'GOBP_VITAMIN_BIOSYNTHETIC_PROCESS',
                                                        'GOBP_VITAMIN_CATABOLIC_PROCESS',
                                                        'GOBP_NAD_METABOLIC_PROCESS',
                                                        'GOBP_NAD_TRANSMEMBRANE_TRANSPORT',
                                                        'GOBP_NAD_TRANSPORT',
                                                        'GOMF_NAD_BINDING',
                                                        'GOBP_RETINAL_METABOLIC_PROCESS',
                                                        'GOBP_CHRONIC_INFLAMMATORY_RESPONSE',
                                                        'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                                                        'GOBP ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS',
                                                        'GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE', 
                                                        'GOBP_CYTOKINE_MEDIATED_SIGNALING_PAHTWAY', 
                                                        'GOBP_CYTOKINE_PRODUCTION', 
                                                        'GOBP_CYTOKINE_PRODUCTION INVOLVED_IN_IMMUNE_RESPONSE',
                                                        'GOBP_CYTOKINE_PRODUCTION_INVLOVED_IN_INFLAMMATORY_RESPONSE',
                                                        'GOBP_RETINAL_CELL_APOPTOTIC_PROCESS',
                                                        'GOBP_RETINAL_GANGLION_CELL_AXON_GUIDANCE',
                                                        'GOBP_ESTABLISHMENT_OF_BLOOD_BARRIER',
                                                        'GOBP_RETINA_HOMEOSTASIS',
                                                        'HP_ABNORMAL_RETINAL_NERVE_FIBER_LAYER_MORPHOLOGY'
                                                        ),species = 'Mus musculus')

ES <- enrichIt(obj = yizhen,
               gene.sets = gene.sets1,
               groups = 1000)

ES2 <- data.frame(yizhen[[]], Idents(yizhen))
colnames(ES2)[ncol(ES2)] <- "cluster"
yizhen <- AddMetaData(yizhen, ES)

ridgeEnrichment(ES2, gene.set = colnames(ES2[12]), group = 'Group', add.rug = TRUE) + facet_wrap(~ident_v2)
