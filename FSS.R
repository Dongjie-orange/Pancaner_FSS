### ---------------
###
### Create: Dongjie Chen
### Date: 2024-04-03 11:25:32
### Email: chen_dj@sjtu.edu.cn
### Pancreatic Disease Center, Ruijin Hospital, SHSMU, Shanghai, China. 
###
### ---------------


rm(list = ls()); gc()
options(stringsAsFactors = F)


#----------------------------------------------------------------------------------
# Step 1: calculate FSx 
#----------------------------------------------------------------------------------

# read the list of cellular senescence-related genes
cellular_senescence_gs <- readRDS("cellular_senescence_gene.rds")

# read the list of pan-cancer scRNA-seq datasets
scRNA_list <- readRDS("40_pancancer_scRNA_seq_dat.rds")

# calculate the FSx
FSx_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    # x <- 1
    sce <- scRNA_list[[x]]
    Idents(sce) <- "Celltype..malignancy."
    sce <- subset(sce, idents = "Malignant cells")
    counts <- sce@assays$RNA@data
    S_CAF_score <- gsva(
      expr = counts, gset.idx.list = cellular_senescence_gs, kcdf = "Gaussian",
      parallel.sz = 60
    )
    S_CAF_score <- as.data.frame(t(getAUC(S_CAF_score)))
    lmGenes <- data.frame(
      gene = rownames(counts),
      coef = NA, p = NA
    )
    for (i in 1:nrow(counts)) {
      cor <- cor.test(counts[i, ], S_CAF_score$cellular_senescence_gs, method = "spearman")
      lmGenes$coef[i] <- cor$estimate
      lmGenes$p[i] <- cor$p.value
    }
    lmGenes$p.adjust <- p.adjust(lmGenes$p, method = "BH")
    return(lmGenes)
  }
)
names(FSx_list) <- names(scRNA_list)





#----------------------------------------------------------------------------------
# Step 2: calculate FSy
#----------------------------------------------------------------------------------

FSy_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    # x <- 1
    sce <- scRNA_list[[x]]
    sce$group <- ifelse(
      sce$Celltype..malignancy. == "Malignant cells", "Malignant",
      "control"
    )
    Idents(sce) <- "group"
    future::plan("multisession", workers = 60)
    DE <- FindMarkers(
      sce, ident.1 = "Malignant", group.by = "group", logfc.threshold = 0.25,
      min.pct = 0.1, base = exp(1)
    )
    DE <- DE[DE$p_val_adj < 1e-05, ]
    return(DE)
  }
)

names(FSy_list) <- names(scRNA_list)





#----------------------------------------------------------------------------------
# Step 3: LM.SIG construction
#----------------------------------------------------------------------------------

identical(
  names(FSx_list),
  names(FSy_list)
)

ls_FSn <- pbapply::pblapply(
  1:length(FSx_list),
  FUN = function(x) {
    # x <- 1
    FSx <- FSx_list[[x]]
    FSy <- FSy_list[[x]]
    FSx <- FSx[FSx$coef > 0 & FSx$p.adjust < 1e-05, ]
    FSy <- rownames(FSy[FSy$avg_logFC > 0.25, ])
    FSy <- FSy[!grepl("^RP[SL]", FSy, ignore.case = F)]  # (ribosome protein free)
    FSn <- FSx[FSx$gene %in% FSy, ]
  }
)

allGenes <- Reduce(rbind, ls_FSn)
allGenes <- unique(allGenes$gene)
allGenesDf <- data.frame(gene = allGenes)

for (i in 1:length(ls_FSn)) {
  allGenesDf <- left_join(
    allGenesDf, ls_FSn[[i]][, c("gene", "coef")],
    by = "gene"
  )
}

allGenesDf <- allGenesDf[!is.na(allGenesDf$gene),
]
rownames(allGenesDf) <- allGenesDf$gene
allGenesDf <- allGenesDf[, -1]
colnames(allGenesDf) <- names(LMx_list)
genelist <- allGenesDf
genelist$all_gmean <- compositions::geometricmeanRow(genelist[, 1:length(ls_FSn)])  # spearmanR geometric mean

sig <- genelist[genelist$all_gmean > 0.25, ]  # filter genes with spearmanR geometric mean > 0.25
sig <- sig[order(sig$all_gmean, decreasing = T),]
FSS <- rownames(sig)

saveRDS(FSS, file = "FSS.rds")




