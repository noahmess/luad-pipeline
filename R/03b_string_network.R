# =============================================================================
# 03b_string_network.R  (Path A)
# Reproducible STRING PPI-network query for the 16 transcriptome-wide DEGs.
# Reads the committed DE table, keeps padj < 0.05, and queries the STRING v12
# REST API (default medium confidence 0.4, Homo sapiens 9606) for:
#   - PPI-enrichment summary (nodes, edges, avg degree, clustering, p-value)
#   - the edge list (to record network composition)
#
# Inputs (committed):
#   outputs/tables/Differential_expression_transcriptomewide.csv
# Outputs (committed):
#   outputs/tables/STRING_network_16DEG.csv     (summary)
#   outputs/tables/STRING_edges_16DEG.csv       (edge list)
#
# NOTE: STRING is a live resource; values reflect the STRING version queried
# (recorded in the summary). Re-running may shift slightly if STRING updates.
# =============================================================================

deg <- read.csv("outputs/tables/Differential_expression_transcriptomewide.csv")
deg16 <- deg$gene[!is.na(deg$padj) & deg$padj < 0.05]
message("03b: querying STRING for ", length(deg16), " DEGs.")

ids     <- paste(deg16, collapse = "%0d")
api     <- "https://string-db.org/api/tsv"
species <- 9606

ppi <- read.delim(url(sprintf("%s/ppi_enrichment?identifiers=%s&species=%d",
                              api, ids, species)))
edges <- read.delim(url(sprintf("%s/network?identifiers=%s&species=%d",
                                api, ids, species)))
ver <- tryCatch(read.delim(url(sprintf("%s/version", api)))$string_version[1],
                error = function(e) NA)

nodes_present <- sort(unique(c(edges$preferredName_A, edges$preferredName_B)))

summary_tbl <- data.frame(
  metric = c("n_DEG_input", "number_of_nodes", "number_of_edges",
             "average_node_degree", "local_clustering_coefficient",
             "expected_number_of_edges", "ppi_enrichment_pvalue",
             "confidence_threshold", "species", "string_version", "source"),
  value = c(length(deg16),
            ppi$number_of_nodes[1],
            ppi$number_of_edges[1],
            ppi$average_node_degree[1],
            ppi$local_clustering_coefficient[1],
            ppi$expected_number_of_edges[1],
            format(ppi$p_value[1], scientific = TRUE),
            "0.4 (medium, STRING default)",
            paste0(species, " (Homo sapiens)"),
            ifelse(is.na(ver), "v12 (unversioned response)", ver),
            "STRING REST API (string-db.org/api/tsv/ppi_enrichment, /network)")
)

write.csv(summary_tbl, "outputs/tables/STRING_network_16DEG.csv", row.names = FALSE)
write.csv(edges[, c("preferredName_A", "preferredName_B", "score")],
          "outputs/tables/STRING_edges_16DEG.csv", row.names = FALSE)

cat("=== STRING PPI-enrichment (16 DEGs) ===\n")
print(summary_tbl, row.names = FALSE)
cat("\nNodes in network (", length(nodes_present), "):\n", sep = "")
cat(paste(nodes_present, collapse = ", "), "\n")
cat("Unique undirected edges:", nrow(edges), "\n")
message("03b: wrote STRING_network_16DEG.csv and STRING_edges_16DEG.csv")
