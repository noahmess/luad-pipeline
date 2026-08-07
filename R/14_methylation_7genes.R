# =============================================================================
# 14_methylation_7genes.R   (version 2)
#
# Figure 14, methylation promotrice des sept genes du panel.
#
# Trois objectifs.
#   1. Regenerer la figure avec LBR, absent de la version actuelle a six panneaux.
#   2. Verifier si le recalcul depuis les donnees 450k concorde avec UALCAN,
#      source annoncee dans les Methodes du manuscrit.
#   3. Diagnostiquer l'origine d'un eventuel ecart, en testant en parallele
#      DEUX jeux de sondes promotrices.
#         jeu A : TSS200 + TSS1500                (jeu UALCAN presume)
#         jeu B : TSS200 + TSS1500 + 5'UTR        (jeu annonce dans tes Methodes)
#      Et DEUX tests, Wilcoxon (script d'origine) et t-test (declare par UALCAN).
#
# La sortie rapporte aussi le nombre de sondes et les effectifs par gene, qui
# peuvent varier selon les valeurs manquantes.
#
# Prerequis. Les fichiers 450k doivent etre en cache sous GDCdata/, ce qui est
# le cas si figure14_methylation_UALCAN.R a deja tourne.
#
# Sorties.
#   outputs/tables/Methylation_7genes_stats.csv
#   outputs/figures/Figure14_methylation_7genes.tiff
# =============================================================================

setwd("~/Desktop/luad-pipeline")

suppressPackageStartupMessages({
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(ggplot2)
})

dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 0. VALEURS UALCAN RELEVEES MANUELLEMENT
# -----------------------------------------------------------------------------
# ualcan.path.uab.edu, TCGA, LUAD, promoter methylation.
# Normal n = 32, tumeur primaire n = 473.
# p exactes, beta lus sur les graphiques donc approximatifs.

ualcan <- data.frame(
    Gene       = c("LMNB1", "LMNB2", "TMPO", "NDC1", "NUP62", "NUP88", "LBR"),
    p_ualcan   = c(5.2086e-01, 1.66861e-01, 1.22550003478494e-09,
                   7.6996e-03, 4.1709e-05, 4.685e-01, 1.1071999983514e-08),
    beta_N_ual = c(0.075, 0.443, 0.085, 0.349, 0.320, 0.062, 0.188),
    beta_T_ual = c(0.069, 0.463, 0.069, 0.366, 0.316, 0.060, 0.174),
    stringsAsFactors = FALSE
)
ualcan$dir_ualcan <- ifelse(ualcan$beta_T_ual > ualcan$beta_N_ual, "hyper", "hypo")

GENES <- ualcan$Gene

# Les deux jeux de sondes a comparer
SETS <- list(
    A_TSSonly = c("TSS1500", "TSS200"),
    B_withUTR = c("TSS1500", "TSS200", "5'UTR", "5UTR")
)


# -----------------------------------------------------------------------------
# 1. CHARGEMENT DES DONNEES 450k
# -----------------------------------------------------------------------------

message("1. Requete GDC, methylation 450k TCGA-LUAD")

query_met <- GDCquery(
    project       = "TCGA-LUAD",
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450"
)

if (!dir.exists("GDCdata/TCGA-LUAD/DNA_Methylation")) {
    message("   cache absent, telechargement, compter 15 a 40 minutes")
    GDCdownload(query_met)
} else {
    message("   cache detecte sous GDCdata/, pas de telechargement")
}

met <- GDCprepare(query_met, summarizedExperiment = TRUE)
message("   objet prepare, ", ncol(met), " echantillons, ", nrow(met), " sondes")


# -----------------------------------------------------------------------------
# 2. ANNOTATION DES SONDES
# -----------------------------------------------------------------------------

rd <- as.data.frame(rowData(met))

col_gene <- grep("Gene_Symbol|gene", colnames(rd), ignore.case = TRUE, value = TRUE)[1]
col_grp  <- grep("Group|Region", colnames(rd), ignore.case = TRUE, value = TRUE)[1]

if (is.na(col_gene) || is.na(col_grp)) {
    cat("\nColonnes disponibles dans rowData :\n"); print(colnames(rd))
    stop("Colonnes d'annotation introuvables. Adapter col_gene et col_grp.")
}
message("2. Annotation, colonne gene = ", col_gene, ", colonne region = ", col_grp)

# Pre-decoupage des annotations, une seule fois, pour la vitesse
gene_split <- strsplit(as.character(rd[[col_gene]]), ";")
grp_split  <- strsplit(as.character(rd[[col_grp]]),  ";")

# Une sonde est retenue pour un gene si au moins une paire (gene, region)
# correspond a la fois au gene cible et a une region du jeu considere.
probes_for <- function(g, regions) {
    keep <- vapply(seq_along(gene_split), function(i) {
        gs <- trimws(gene_split[[i]]); rs <- trimws(grp_split[[i]])
        if (length(gs) == 0 || length(rs) == 0) return(FALSE)
        n <- max(length(gs), length(rs))
        gs <- rep(gs, length.out = n); rs <- rep(rs, length.out = n)
        any(gs == g & rs %in% regions)
    }, logical(1))
    rownames(rd)[keep]
}

probe_map <- list()
for (sn in names(SETS)) {
    probe_map[[sn]] <- lapply(GENES, function(g) probes_for(g, SETS[[sn]]))
    names(probe_map[[sn]]) <- GENES
}

cat("\n=== Nombre de sondes promotrices par gene ===\n")
cnt <- data.frame(
    Gene      = GENES,
    n_probes_A = sapply(GENES, function(g) length(probe_map$A_TSSonly[[g]])),
    n_probes_B = sapply(GENES, function(g) length(probe_map$B_withUTR[[g]])),
    stringsAsFactors = FALSE
)
cnt$diff_5UTR <- cnt$n_probes_B - cnt$n_probes_A
print(cnt, row.names = FALSE)
cat("\nA = TSS200 + TSS1500, jeu UALCAN presume\n")
cat("B = A + 5'UTR, jeu annonce dans les Methodes actuelles\n")
cat("diff_5UTR = sondes ajoutees par l'inclusion du 5'UTR\n")
if (all(cnt$diff_5UTR == 0)) {
    cat("\nAucune sonde 5'UTR supplementaire. Les deux jeux sont identiques,\n")
    cat("la question de la definition des sondes est donc sans objet ici.\n")
} else {
    cat("\nLes deux jeux different. Comparer les p-values des deux jeux\n")
    cat("ci-dessous pour savoir si l'ecart avec UALCAN vient de la.\n")
}


# -----------------------------------------------------------------------------
# 3. MOYENNE PAR GENE ET TESTS, POUR CHAQUE JEU DE SONDES
# -----------------------------------------------------------------------------

beta <- assay(met)
grp_sample <- ifelse(substr(colnames(beta), 14, 15) == "11", "Normal", "Tumor")

cat("\n=== Effectifs globaux ===\n")
cat("Normal", sum(grp_sample == "Normal"), "  Tumeur", sum(grp_sample == "Tumor"), "\n")
cat("UALCAN annonce Normal 32 et Tumeur 473.\n")

stats_for_set <- function(sn) {
    pm <- probe_map[[sn]]
    do.call(rbind, lapply(GENES, function(g) {
        p <- pm[[g]]
        if (length(p) == 0) {
            return(data.frame(Gene = g, set = sn, n_probes = 0,
                              n_N = 0, n_T = 0, beta_N = NA, beta_T = NA,
                              delta = NA, p_wilcox = NA, p_ttest = NA,
                              direction = NA, stringsAsFactors = FALSE))
        }
        v <- if (length(p) == 1) as.numeric(beta[p, ]) else colMeans(beta[p, , drop = FALSE], na.rm = TRUE)
        n  <- v[grp_sample == "Normal"]; n  <- n[!is.na(n)]
        t_ <- v[grp_sample == "Tumor"];  t_ <- t_[!is.na(t_)]
        if (length(n) < 3 || length(t_) < 3) {
            return(data.frame(Gene = g, set = sn, n_probes = length(p),
                              n_N = length(n), n_T = length(t_),
                              beta_N = NA, beta_T = NA, delta = NA,
                              p_wilcox = NA, p_ttest = NA, direction = NA,
                              stringsAsFactors = FALSE))
        }
        mN <- median(n); mT <- median(t_)
        data.frame(
            Gene      = g,
            set       = sn,
            n_probes  = length(p),
            n_N       = length(n),
            n_T       = length(t_),
            beta_N    = round(mN, 4),
            beta_T    = round(mT, 4),
            delta     = round(mT - mN, 4),
            p_wilcox  = signif(wilcox.test(t_, n)$p.value, 4),
            p_ttest   = signif(t.test(t_, n, var.equal = FALSE)$p.value, 4),
            direction = ifelse(mT > mN, "hyper", "hypo"),
            stringsAsFactors = FALSE
        )
    }))
}

resA <- stats_for_set("A_TSSonly")
resB <- stats_for_set("B_withUTR")

cat("\n=== Jeu A, TSS200 + TSS1500 ===\n")
print(resA[, c("Gene","n_probes","n_N","n_T","beta_N","beta_T","delta",
               "p_wilcox","p_ttest","direction")], row.names = FALSE)

cat("\n=== Jeu B, avec 5'UTR ===\n")
print(resB[, c("Gene","n_probes","n_N","n_T","beta_N","beta_T","delta",
               "p_wilcox","p_ttest","direction")], row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. CONCORDANCE AVEC UALCAN, POUR CHAQUE JEU
# -----------------------------------------------------------------------------

sigf <- function(p) ifelse(is.na(p), NA, p < 0.05)

concord <- function(res, label) {
    m <- merge(res, ualcan, by = "Gene", sort = FALSE)
    m <- m[match(GENES, m$Gene), ]
    m$dir_ok <- ifelse(m$direction == m$dir_ualcan, "oui", "NON")
    m$sig_ok <- ifelse(sigf(m$p_ttest) == sigf(m$p_ualcan), "oui", "NON")
    cat("\n=== Concordance UALCAN, jeu", label, "===\n")
    print(m[, c("Gene","n_probes","delta","direction","p_ttest",
                "p_ualcan","dir_ok","sig_ok")], row.names = FALSE)
    c(dir = sum(m$dir_ok == "NON", na.rm = TRUE),
      sig = sum(m$sig_ok == "NON", na.rm = TRUE))
}

cA <- concord(resA, "A (TSS only)")
cB <- concord(resB, "B (avec 5'UTR)")

cat("\n=== VERDICT ===\n")
cat("Jeu A, discordances :", cA["dir"], "direction,", cA["sig"], "significativite\n")
cat("Jeu B, discordances :", cB["dir"], "direction,", cB["sig"], "significativite\n\n")

best <- if (sum(cA) <= sum(cB)) "A" else "B"
if (sum(cA) == 0 || sum(cB) == 0) {
    cat("Au moins un jeu reproduit UALCAN. Retenir le jeu", best,
        "et aligner la phrase des Methodes sur sa definition.\n")
} else if (cA["dir"] == 0 && cB["dir"] == 0) {
    cat("Les directions concordent dans les deux jeux, seules des\n")
    cat("significativites divergent. Retenir le jeu", best, "et afficher les\n")
    cat("p-values UALCAN dans la legende, en precisant que les distributions\n")
    cat("sont recalculees.\n")
} else {
    cat("Divergence de direction detectee. Ne pas ecrire le texte avant\n")
    cat("d'avoir identifie la cause. Verifier l'annotation des sondes.\n")
}

# Amplitude maximale parmi les genes significatifs selon UALCAN
resBest <- if (best == "A") resA else resB
mm <- merge(resBest, ualcan, by = "Gene", sort = FALSE)
sig_genes <- mm[mm$p_ualcan < 0.05 & !is.na(mm$delta), ]
if (nrow(sig_genes)) {
    cat("\nAmplitude, genes significatifs selon UALCAN, jeu", best, "\n")
    print(sig_genes[, c("Gene","delta","direction","p_ualcan")], row.names = FALSE)
    cat("Delta absolu maximal :", max(abs(sig_genes$delta)), "\n")
    cat("Ce chiffre alimente l'argument d'amplitude de la section 3.7.\n")
    nd1 <- sig_genes[sig_genes$Gene == "NDC1", ]
    if (nrow(nd1)) {
        cat("NDC1, delta =", nd1$delta, ", direction", nd1$direction,
            ". A confirmer positif avant d'ecrire la discordance directionnelle.\n")
    }
}

out <- rbind(resA, resB)
write.csv(out, "outputs/tables/Methylation_7genes_stats.csv", row.names = FALSE)
cat("\nEcrit : outputs/tables/Methylation_7genes_stats.csv\n")


# -----------------------------------------------------------------------------
# 5. FIGURE A SEPT PANNEAUX
# -----------------------------------------------------------------------------
# Construite sur le jeu retenu. Panneaux ordonnes par p-value UALCAN croissante.

pm_best <- probe_map[[if (best == "A") "A_TSSonly" else "B_withUTR"]]

mat <- t(sapply(GENES, function(g) {
    p <- pm_best[[g]]
    if (length(p) == 0) return(rep(NA_real_, ncol(beta)))
    if (length(p) == 1) return(as.numeric(beta[p, ]))
    colMeans(beta[p, , drop = FALSE], na.rm = TRUE)
}))
colnames(mat) <- colnames(beta)

ordre <- ualcan$Gene[order(ualcan$p_ualcan)]
genes_ok <- ordre[ordre %in% resBest$Gene[!is.na(resBest$beta_N)]]

df <- do.call(rbind, lapply(genes_ok, function(g) {
    data.frame(Gene = g, Group = grp_sample, Beta = as.numeric(mat[g, ]),
               stringsAsFactors = FALSE)
}))
df <- df[!is.na(df$Beta), ]
df$Gene  <- factor(df$Gene, levels = genes_ok)
df$Group <- factor(df$Group, levels = c("Normal", "Tumor"))

lab <- setNames(
    sprintf("%s\n(p = %s)", ualcan$Gene,
            format(ualcan$p_ualcan, digits = 3, scientific = TRUE)),
    ualcan$Gene
)

p <- ggplot(df, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 2,
               labeller = labeller(Gene = lab)) +
    scale_fill_manual(values = c(Normal = "#4C72B0", Tumor = "#C44E52")) +
    labs(x = NULL, y = "Promoter methylation (beta value)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank())

ggsave("outputs/figures/Figure14_methylation_7genes.tiff", p,
       width = 11, height = 6.5, dpi = 300, compression = "lzw")

cat("Ecrit : outputs/figures/Figure14_methylation_7genes.tiff\n")
cat("Jeu de sondes utilise pour la figure :", best, "\n")
cat("\nLegende a prevoir, transparence complete sur les deux sources.\n")
cat("  Beta values recomputed from TCGA Illumina 450k data.\n")
cat("  p-values from UALCAN.\n")

cat("\n=== Fin ===\n")
