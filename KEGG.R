# ============================================================
# KEGG pathway analysis (GSEA + ORA), HIV vs Mock, all clusters
# Uses the stored per-cluster DE files (positive avg_log2FC = higher in HIV).
# Ribosomal and mitochondrial genes are removed before the analysis.
# Needs an internet connection once (downloads the KEGG human pathway sets).
# ============================================================
library(dplyr); library(tidyr); library(ggplot2)
library(clusterProfiler); library(org.Hs.eg.db)

# ---- Settings -------------------------------------------------
out_root <- "E:/ssDNA_vs_Mock_GO"
kegg_dir <- file.path(out_root, "KEGG")

obj_path      <- "E:/multiome_with_Marie_celltype_atac.rds"
rna_assay     <- "RNA"
de_dir   <- file.path(out_root, "DE")
kegg_dir <- file.path(out_root, "KEGG")
cluster_col   <- "seurat_clusters"
condition_col <- "condition"
hiv_label     <- "ssDNA"
mock_label    <- "Mock"
padj_cut <- 0.05; lfc_cut <- 0.25
n_show   <- 10
fs       <- 16
dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)

# Genes removed before GSEA and ORA (the saved DE files are NOT changed):
#   MT-*            mitochondrial-encoded genes
#   RPL*/RPS*       cytosolic ribosomal proteins (incl. RPLP0/1/2, RPSA, ribosomal pseudogenes)
#   MRPL*/MRPS*     mitochondrial ribosomal proteins (remove this part if you want to keep them)
drop_pattern <- paste0(
  "^MT-|",                  # Mitochondrial genes
  "^RP[SL][0-9]|",           # RPS/RPL ribosomal genes
  "^RPLP[0-2]$|",            # RPLP0, RPLP1, RPLP2
  "^RPSA$|",                 # RPSA
  "^MRP[SL][0-9]|",          # Mitochondrial ribosomal genes
  "^MALAT1$|",               # MALAT1
  "^(XIST|TSIX)$|",          # X-linked genes
  "^(UTY|USP9Y|PRKY|DDX3Y|EIF1AY)$|", # Y-linked genes
  "^RPS4Y[12]$|",            # RPS4Y1, RPS4Y2
  "^ZFY$|",                  # Y-linked gene
  "^KDM5D$|",                # Y-linked gene
  "^DDX3X$|",                # X-linked gene
  "^EIF1AX$|",               # X-linked gene
  "^RPS4X$|",                # X-linked gene
  "^KDM6A$|",                # X-linked gene
  "^UTX$"                    # Alternative symbol for KDM6A
)

lab_hiv  <- paste("Higher in", hiv_label)
lab_mock <- paste("Higher in", mock_label)
dir_cols <- setNames(c("firebrick", "steelblue"), c(lab_hiv, lab_mock))

# ---- KEGG gene sets (downloaded once) --------------------------
message("Downloading KEGG human pathways...")
kegg <- download_KEGG("hsa")
t2g <- data.frame(term = kegg$KEGGPATHID2EXTID$from, gene = kegg$KEGGPATHID2EXTID$to)
t2n <- data.frame(term = kegg$KEGGPATHID2NAME$from,
                  name = sub(" - Homo sapiens \\(human\\)$", "", kegg$KEGGPATHID2NAME$to))
write.csv(t2n, file.path(kegg_dir, "KEGG_hsa_pathway_names.csv"), row.names = FALSE)

to_entrez <- function(symbols) {
  m <- suppressMessages(bitr(unique(symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  m[!duplicated(m$SYMBOL) & !duplicated(m$ENTREZID), ]
}

bar_plot <- function(d, xcol, title, xlab, fn, w = 10) {
  d$Description <- factor(d$Description, levels = unique(d$Description[order(d[[xcol]])]))
  p <- ggplot(d, aes(.data[[xcol]], Description, fill = direction)) + geom_col() +
    scale_fill_manual(values = dir_cols) +
    scale_y_discrete(labels = scales::label_wrap(45)) +
    labs(title = title, x = xlab, y = NULL, fill = NULL) +
    theme_classic(base_size = fs) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(fn, p, width = w, height = max(4, 0.35 * nrow(d) + 2), dpi = 300, bg = "white")
}

# ---- Per cluster ------------------------------------------------
files <- list.files(de_dir, pattern = "^Cluster_.*\\.csv$")
cat("DE files:", length(files), "\n")
if (length(files) == 0) stop("No DE files found in ", de_dir)
gsea_all <- list(); ora_all <- list()

for (f in files) {
  cl <- sub(paste0("_HIV_vs_", mock_label, "$"), "", sub("^Cluster_", "", sub("\\.csv$", "", f)))
  message("\nCluster ", cl)
  deg <- read.csv(file.path(de_dir, f), check.names = FALSE)
  cdir <- file.path(kegg_dir, paste0("Cluster_", cl)); dir.create(cdir, showWarnings = FALSE)
  
  # remove ribosomal / mitochondrial genes
  n0  <- nrow(deg)
  deg <- deg[!grepl(drop_pattern, deg$gene, ignore.case = TRUE), ]
  message("Removed ", n0 - nrow(deg), " ribosomal/mitochondrial genes; ", nrow(deg), " remain")
  
  map <- to_entrez(deg$gene)
  deg <- deg %>% inner_join(map, by = c("gene" = "SYMBOL")) %>% filter(!is.na(avg_log2FC))
  
  # ---- GSEA ----
  rk <- sort(setNames(deg$avg_log2FC, deg$ENTREZID), decreasing = TRUE)
  g <- tryCatch(
    GSEA(rk, TERM2GENE = t2g, TERM2NAME = t2n, minGSSize = 10, maxGSSize = 500,
         pvalueCutoff = 1, eps = 0, verbose = FALSE),
    error = function(e) { message("GSEA failed: ", e$message); NULL })
  
  if (!is.null(g) && nrow(g@result) > 0) {
    g <- tryCatch(setReadable(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(e) g)
    res <- g@result %>% mutate(cluster = cl) %>% arrange(desc(NES))
    write.csv(res, file.path(cdir, "KEGG_GSEA.csv"), row.names = FALSE)
    gsea_all[[cl]] <- res
    
    sig <- res %>% filter(p.adjust < padj_cut)
    if (nrow(sig) > 0) {
      top <- bind_rows(sig %>% filter(NES > 0) %>% slice_max(NES, n = n_show, with_ties = FALSE),
                       sig %>% filter(NES < 0) %>% slice_min(NES, n = n_show, with_ties = FALSE)) %>%
        mutate(direction = ifelse(NES > 0, lab_hiv, lab_mock))
      bar_plot(top, "NES", paste0("Cluster ", cl, ": KEGG GSEA"), "NES",
               file.path(cdir, "KEGG_GSEA_bar.png"))
    }
  }
  
  # ---- ORA (up / down DEGs) ----
  sigdeg <- deg %>% filter(p_val_adj < padj_cut, abs(avg_log2FC) >= lfc_cut) %>%
    mutate(direction = ifelse(avg_log2FC > 0, lab_hiv, lab_mock))
  ora_cl <- list()
  for (d in c(lab_hiv, lab_mock)) {
    genes <- sigdeg$ENTREZID[sigdeg$direction == d]
    if (length(genes) < 5) next
    e <- tryCatch(
      enricher(genes, TERM2GENE = t2g, TERM2NAME = t2n, universe = deg$ENTREZID,
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10),
      error = function(e) { message("ORA failed (", d, "): ", e$message); NULL })
    if (is.null(e) || nrow(e@result[e@result$p.adjust < 0.05, ]) == 0) next
    e <- tryCatch(setReadable(e, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), error = function(err) e)
    ora_cl[[d]] <- e@result %>% filter(p.adjust < 0.05) %>%
      mutate(direction = d, cluster = cl) %>% arrange(p.adjust)
  }
  if (length(ora_cl) > 0) {
    ora <- bind_rows(ora_cl)
    write.csv(ora, file.path(cdir, "KEGG_ORA.csv"), row.names = FALSE)
    ora_all[[cl]] <- ora
    top <- ora %>% group_by(direction) %>% slice_min(p.adjust, n = 15, with_ties = FALSE) %>% ungroup() %>%
      mutate(neglog10 = -log10(p.adjust))
    bar_plot(top, "neglog10", paste0("Cluster ", cl, ": KEGG ORA"), expression(-log[10]("adj. P")),
             file.path(cdir, "KEGG_ORA_bar.png"))
  }
}

# ---- Combined tables + heatmap ----------------------------------
if (length(gsea_all) == 0) stop("No KEGG GSEA results; check the DE folder and internet connection.")
G <- bind_rows(gsea_all)
write.csv(G, file.path(kegg_dir, "ALL_KEGG_GSEA.csv"), row.names = FALSE)
write.csv(G %>% filter(p.adjust < padj_cut) %>% arrange(cluster, desc(NES)),
          file.path(kegg_dir, "ALL_KEGG_GSEA_significant.csv"), row.names = FALSE)
if (length(ora_all) > 0)
  write.csv(bind_rows(ora_all), file.path(kegg_dir, "ALL_KEGG_ORA.csv"), row.names = FALSE)

# heatmap: pathways significant in at least one cluster (top 40 by max |NES|)
keep <- G %>% filter(p.adjust < padj_cut) %>% group_by(Description) %>%
  summarise(m = max(abs(NES))) %>% slice_max(m, n = 40, with_ties = FALSE) %>% pull(Description)
if (length(keep) > 0) {
  H <- G %>% filter(Description %in% keep)
  cl_levels <- unique(G$cluster)
  cl_levels <- cl_levels[order(suppressWarnings(as.numeric(cl_levels)), cl_levels)]
  H <- H %>% mutate(cluster = factor(cluster, levels = cl_levels),
                    star = ifelse(p.adjust < padj_cut, "*", ""))
  ord <- H %>% group_by(Description) %>% summarise(m = mean(NES, na.rm = TRUE)) %>% arrange(m) %>% pull(Description)
  H$Description <- factor(H$Description, levels = ord)
  
  p <- ggplot(H, aes(cluster, Description, fill = NES)) + geom_tile(colour = "white") +
    geom_text(aes(label = star), size = 6, vjust = 0.75) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
    labs(title = paste("KEGG GSEA:", hiv_label, "vs", mock_label),
         subtitle = paste0("Red = higher in ", hiv_label, "; * = adj. P < ", padj_cut),
         x = "Cluster", y = NULL) +
    theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(kegg_dir, "KEGG_NES_heatmap.png"), p,
         width = max(8, 0.5 * length(cl_levels) + 7), height = max(6, 0.3 * length(keep) + 3),
         dpi = 300, bg = "white")
  ggsave(file.path(kegg_dir, "KEGG_NES_heatmap.pdf"), p,
         width = max(8, 0.5 * length(cl_levels) + 7), height = max(6, 0.3 * length(keep) + 3))
}

message("Done. Results in: ", kegg_dir)

