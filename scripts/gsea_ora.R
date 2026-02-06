###################### GSEA + ORA: NM Protein-coding ######################

# 1. SETUP
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DESeq2, clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, 
               dplyr, msigdbr, ggplot2, tibble, pathview, ggupset)

# 2. DESeq2
path_counts <- "/Users/magali/server/tmp/Pruebas/rebe/bowtie2/bam/BX_count_matrix.txt"
counts <- read.delim(path_counts, row.names=1)
design <- data.frame(condition = factor(c("AZ","AZ","AZ","Control","Control","Control"), 
                                        levels=c("Control","AZ")), row.names=colnames(counts))

dds <- DESeqDataSetFromMatrix(counts, design, ~condition) |> DESeq()
deseq_stats <- results(dds)

# 3. NM protein-coding + anotación
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

cat(sprintf("NM: %d total | ENTREZ: %d\n", nrow(nm_proteins), sum(!is.na(nm_proteins$entrez_id))))

# 4. 1 transcrito/gen (mayor |stat|, desempate |log2FC|)
gsea_candidates <- nm_proteins |> 
  filter(!is.na(entrez_id)) |> 
  slice_max(abs(stat), by=entrez_id, with_ties=FALSE) |> 
  slice_max(abs(log2FoldChange), by=entrez_id, with_ties=FALSE)

cat(sprintf("GSEA genes únicos: %d\n", nrow(gsea_candidates)))

# 5. Ranking GSEA
gsea_ranking <- gsea_candidates$stat
names(gsea_ranking) <- make.unique(gsea_candidates$entrez_id)
gsea_ranking <- sort(gsea_ranking[!is.na(names(gsea_ranking))], decreasing=TRUE)
#head(gsea_ranking, 10)

# 6. ORA: Genes DE
de_up_table <- nm_proteins |>
  filter(padj < 0.05 & log2FoldChange > 1) |>
#  unique(nm_proteins$entrez_id[nm_proteins$padj < 0.05 & nm_proteins$log2FoldChange > 1]) |> 
  mutate(
    entrez_id = mapIds(org.Hs.eg.db, refseq_base, "ENTREZID", "REFSEQ"),
    gene_sym  = mapIds(org.Hs.eg.db, refseq_base, "SYMBOL", "REFSEQ")) |>
  filter (!is.na(entrez_id))

de_up_unique <- de_up_table |> filter(padj < 0.05,
                                      log2FoldChange > 1,
                                      !is.na(gene_sym),
                                      !is.na(entrez_id)
                                      ) 
de_up_unique <-  unique(de_up_unique$entrez_id)
de_up_unique <- as.data.frame(de_up_unique, head =TRUE, sep ="\t")


de_down <- unique(nm_proteins$entrez_id[nm_proteins$padj < 0.05 & nm_proteins$log2FoldChange < -1])
gene_universe <- unique(nm_proteins$entrez_id)

# 7. Análisis funcionales
set.seed(1234)

# ORA KEGG UP
kegg_ora_up <- enrichKEGG(de_up, universe=gene_universe, organism="hsa", pvalueCutoff =0.05, qvalueCutoff = 0.05) |> 
  setReadable (OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
                                
kegg_ora_up_table <- as.data.frame(kegg_ora_up)

upsetplot(kegg_ora_up)

heatplot(kegg_ora_up, showCategory = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

kegg_sim <- pairwise_termsim(kegg_ora_up)

treeplot(kegg_sim, 
         showCategory = 15,
         nCluster = 3)

# Para crear el mapa de la vía

# 1. Crear el vector
pathview_input <- nm_proteins$log2FoldChange
names(pathview_input) <- nm_proteins$entrez_id

# 2. Aplicar el filtro usando los nombres para evitar desajustes
significativos <- nm_proteins$entrez_id[nm_proteins$padj < 0.05 & !is.na(nm_proteins$padj)]
pathview_input[!(names(pathview_input) %in% significativos)] <- 0

# 3. Limpiar NAs de los nombres
pathview_input <- pathview_input[!is.na(names(pathview_input))]

pathview(gene.data  = pathview_input,
         pathway.id = "hsa04140", # ID de la vía (ej. Viral protein interaction)
         species    = "hsa",
         kegg.dir   = ".",        # Directorio donde se guarda
         low = list(gene = "blue"), # Color para genes que bajan
         high = list(gene = "red"), # Color para genes que suben
         mid = list(gene = "white")) # Genes que se detectaron pero no mantienen cambios



#kegg_ora_up_tabla <- as.data.frame(kegg_ora_up, header = TRUE, sep = "\t") |>


#kegg_ora_up@result |> slice_head(n=10)

dotplot(kegg_ora_up, 
        showCategory = 10, 
        color = "qvalue",
          orderBy = "GeneRatio") + 
  scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "q-value") +
  labs(title = "Enriched KEGG Pathways", 
       x = "Gene Ratio", 
       y = "") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank())

dotplot(kegg_ora_up, 
        showCategory = 10, 
        color = "qvalue",
          orderBy = "GeneRatio") + 
  scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "q-value") +
  labs(title = "Enriched KEGG Pathways", 
       x = "Gene Ratio", 
       y = "") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank())



# GSEA KEGG
kegg_gsea <- gseKEGG(gsea_ranking, organism="hsa", minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)
dotplot(kegg_gsea, 
        showCategory = 10, 
        split = ".sign",           # Movido dentro de dotplot
        color = "p.adjust",       
        orderBy = "GeneRatio") +   
  facet_grid(. ~ .sign) +          # Organiza por Activados/Reprimidos
  scale_color_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adj") +
  labs(title = "GSEA KEGG Pathways", 
       x = "Gene Ratio", 
       y = "") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

kegg_df <- as.data.frame(kegg_gsea) |>
  filter(p.adjust < 0.05) |>
  group_by(sign = ifelse(NES > 0, "Activated", "Suppressed")) %>%
  slice_max(abs(NES), n = 10) |>
  ungroup()

ggplot(kegg_df, aes(x = NES, y = reorder(Description, NES), fill = p.adjust)) +
  geom_col() +
  # Dividir en paneles por Activado/Reprimido
  facet_grid(sign ~ ., scales = "free_y") +
  scale_fill_continuous(low = "#E41A1C", high = "#377EB8", name = "p.adj") +
  labs(title = "GSEA KEGG Pathways",
       x = "Normalized Enrichment Score (NES)", 
       y = NULL) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# GSEA GO BP
go_gsea <- gseGO(gsea_ranking, org.Hs.eg.db, ont="BP", minGSSize=10, pvalueCutoff=0.05)
dotplot(go_gsea, showCategory=10,split=".sign") + facet_grid(.~.sign)

# GSEA Hallmarks
hallmark_sets <- msigdbr(species="human", category="H")[, c("gs_name", "entrez_gene")]
colnames(hallmark_sets) <- c("pathway", "entrez_gene")
hallmark_gsea <- GSEA(gsea_ranking, TERM2GENE=hallmark_sets, minGSSize=10, pvalueCutoff=0.05)
dotplot(hallmark_gsea, showCategory=10, split=".sign") + facet_grid(.~.sign)

# 8. Pathview (vía específica)
pathview(gene.data=gsea_ranking, pathway.id="hsa04110", species="hsa")

# 9. Tablas legibles
kegg_gsea_readable <- setReadable(kegg_gsea, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
go_gsea_readable <- setReadable(go_gsea, OrgDb=org.Hs.eg.db, keyType="ENTREZID")

# 10. Exportar
write.csv(gsea_candidates, "01_gsea_input.csv", row.names=FALSE)
write.csv(kegg_gsea_readable@result, "02_gsea_kegg.csv", row.names=FALSE)
write.csv(go_gsea_readable@result, "03_gsea_go.csv", row.names=FALSE)
write.csv(hallmark_gsea@result, "04_gsea_hallmarks.csv", row.names=FALSE)
write.csv(kegg_ora_up@result, "05_ora_kegg.csv", row.names=FALSE)
