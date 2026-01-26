###################### GSEA + ORA (NM Protein-coding) ######################

# ===== 0. SETUP =====
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DESeq2, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, 
               dplyr, msigdbr, ggplot2, tibble)

cat("✓ Setup completo\n")

# ===== 1. DESeq2: Counts → Estadísticos =====
path_counts <- "/Users/magali/server/tmp/Pruebas/rebe/bowtie2/bam/BX_count_matrix.txt"
counts <- read.delim(path_counts, row.names=1)
design <- data.frame(condition = factor(c("AZ","AZ","AZ","Control","Control","Control"), 
                                        levels=c("Control","AZ")), 
                     row.names=colnames(counts))

dds <- DESeqDataSetFromMatrix(counts, design, ~condition) |> DESeq()
deseq_stats <- results(dds)

cat("✓ DESeq2 completado\n")

# ===== 2. Protein-coding NM + Anotación =====
nm_proteins <- tibble(refseq_id = rownames(deseq_stats)) |> 
  bind_cols(as.data.frame(deseq_stats)) |> 
  filter(grepl("^NM_", refseq_id)) |> 
  mutate(refseq_base = sub("\\..*", "", refseq_id)) |> 
  mutate(
    entrez_id = mapIds(org.Hs.eg.db, refseq_base, "ENTREZID", "REFSEQ"),
    gene_sym  = mapIds(org.Hs.eg.db, refseq_base, "SYMBOL", "REFSEQ")
  ) |> 
  filter(!is.na(stat)) |> 
  arrange(desc(abs(stat)))

cat(sprintf("NM proteins: %d total | %d con ENTREZ\n", 
            nrow(nm_proteins), sum(!is.na(nm_proteins$entrez_id))))

# ===== 3. SELECCIÓN: 1 transcripto por gen =====
# Criterio: Mayor |stat| y desempate por mayor |log2FC|
gsea_candidates <- nm_proteins |> 
  filter(!is.na(entrez_id)) |> 
  slice_max(abs(stat), by=entrez_id, with_ties=FALSE) |> 
  slice_max(abs(log2FoldChange), by=entrez_id, with_ties=FALSE)

cat(sprintf("Genes únicos seleccionados: %d\n", nrow(gsea_candidates)))

# ===== 4. GSEA =====
# Preparar Ranking
gsea_ranking <- gsea_candidates$stat
names(gsea_ranking) <- gsea_candidates$entrez_id
gsea_ranking <- sort(gsea_ranking, decreasing=TRUE)


# Enriquecimiento funcional con KEGG + GSEA (vías apagadas/encendidas globales)

gsea_kegg <- gseKEGG(geneList     = gsea_ranking,
                     organism     = "hsa",
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)

# Visualizar
dotplot(gsea_kegg, showCategory=15, split=".sign") + 
  facet_grid(.~.sign) +
  ggtitle("GSEA KEGG: AZ vs Control")

# Ver la tabla de resultados (convertir IDs a nombres)
gsea_kegg_readable <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
View(as.data.frame(gsea_kegg_readable))


if (!require("pacman")) install.packages("pacman")
pacman::p_load(pathview)

# Esto generará un archivo .png en tu carpeta de trabajo
pathview(gene.data = gsea_ranking, 
         pathway.id = "hsa04110", # El ID de la vía que te interesa
         species    = "hsa", 
         limit      = list(gene=max(abs(gsea_ranking))), 
         low = "blue", high = "red") # Azul = baja, Rojo = sube



# 4a. GSEA GO BP
set.seed(1234)
gsea_go <- gseGO(gsea_ranking, org.Hs.eg.db, ont="BP", keyType="ENTREZID",
                 minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, eps=0)

tabla_gsea <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
tabla_gsea <- as.data.frame(tabla_gsea)
View(tabla_gsea)

busqueda <- tabla_gsea |> filter(grepl("GO:0097305", Description, ignore.case = TRUE))

print(dotplot(gsea_go, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
        ggtitle("GSEA GO BP: AZ vs Control"))

# 4b. GSEA Hallmarks
hallmark_sets <- msigdbr(species="human", category="H") %>% 
  dplyr::select(gs_name, entrez_gene) %>% 
  rename(pathway = gs_name)

gsea_hallmarks <- GSEA(gsea_ranking, TERM2GENE=hallmark_sets,
                       minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)

print(dotplot(gsea_hallmarks, showCategory=15, split=".sign") + facet_grid(.~.sign))

# ===== 5. ORA (Over-Representation Analysis) =====
# Definir listas de genes
de_up <- nm_proteins |> filter(padj < 0.05, log2FoldChange > 1) |> pull(entrez_id) |> unique()

de_down <- nm_proteins |> filter(padj < 0.05, log2FoldChange < -1) |> pull(entrez_id) |> unique()
gene_universe <- unique(nm_proteins$entrez_id)

cat(sprintf("ORA Setup: %d up | %d down | fondo %d\n", length(de_up), length(de_down), length(gene_universe)))

# 5a. ORA GO BP
ora_go_up <- enrichGO(de_up, universe=gene_universe, OrgDb=org.Hs.eg.db, 
                      ont="BP", keyType="ENTREZID", pvalueCutoff=0.05, readable=TRUE)

tabla_ora_go_up <- as.data.frame(ora_go_up)


ora_go_down <- enrichGO(de_down, universe=gene_universe, OrgDb=org.Hs.eg.db, 
                        ont="BP", keyType="ENTREZID", pvalueCutoff=0.05, readable=TRUE)

tabla_ora_go_down <- as.data.frame(ora_go_down)

print(dotplot(ora_go_up, showCategory=15) + ggtitle("ORA GO UP"))
print(dotplot(ora_go_down, showCategory=15) + ggtitle("ORA GO DOWN"))

# ===== ORA KEGG (Genes UP) =====
ora_kegg_up <- enrichKEGG(gene         = de_up,
                          organism     = "hsa",
                          universe     = gene_universe,
                          pvalueCutoff = 0.05)

# Visualizar UP
dotplot(ora_kegg_up, showCategory=10) + ggtitle("ORA KEGG: Genes UP")

# ===== ORA KEGG (Genes DOWN) =====
ora_kegg_down <- enrichKEGG(gene         = de_down,
                            organism     = "hsa",
                            universe     = gene_universe,
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 0.5)


# Visualizar DOWN
dotplot(ora_kegg_down, showCategory=10) + ggtitle("ORA KEGG: Genes DOWN")




# ===== 6. EXPORTAR =====
# Descomentar para guardar archivos:
# write.csv(gsea_candidates, "01_gsea_input_NM.csv", row.names=FALSE)
# write.csv(gsea_go@result, "02_gsea_go_bp.csv", row.names=FALSE)
# write.csv(gsea_hallmarks@result, "03_gsea_hallmarks.csv", row.names=FALSE)
# write.csv(ora_go_up@result, "04_ora_go_up.csv", row.names=FALSE)
# write.csv(ora_kegg_up@result, "05_ora_kegg_up.csv", row.names=FALSE)

