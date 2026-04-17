# ============================================================================
# Script:  Top_intron_match_analysis_NSC.R
# Author: Yali Zhang
# Date: 2026-04-17
#
# Description:
#   Identifies de novo variant (DNV) burden in top 55 TIR introns compared 
#   with size- and GC-matched background introns, using two ASD cohorts
#   (SPARK and SSC).
#
# Workflow:
#   1. Derive intron coordinates for NSC-expressed genes (TPM >= 1 in both
#      replicates) from GENCODE v38 annotation.
#   2. Load top TIR introns and build a background intron set by excluding
#      any introns that overlap top TIR introns.
#   3. For each of 100 iterations, randomly match each TIR intron to a background 
#      intron of similar length (±10 %) and GC content (±5 %).
#   4. Compute per-individual DNV counts in the target vs. matched introns and
#      test group differences with a GEE model accounting for family structure.
#   5. Visualize the distribution of 100 DNV rates and summarise the
#      proportion of iterations reaching each significance threshold.
#
# Input files:
#   - gencode.v38.annotation.gtf
#   - RNA-Seq_iPSC.NSC.Neuron_tpm.txt
#   - top_TIR_Neuron_introns_hg38.txt
#   - SPARK_DNV.tsv
#   - SSC_DNV.tsv
#
# Output:
#   - SPARK_TOP_matched_intron  (ggplot histogram)
#   - SSC_TOP_matched_intron    (ggplot histogram)
#   - sig_fig                   (ggplot bar chart of significance levels)
# ============================================================================

library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(drgee)

genome <- BSgenome.Hsapiens.UCSC.hg38

# ── 1. Prepare genomic intron list ───────────────────────────────────────────

# Load GENCODE v38 GTF annotation and NSC RNA-seq TPM table
gtf <- import("gencode.v38.annotation.gtf")
tpm_table <- read.delim("RNA-Seq_iPSC.NSC.Neuron_tpm.txt", header = TRUE)

# Keep genes with TPM >= 1 in both NSC replicates
expressed_NSC_genes <- tpm_table |>
  filter(NSC_r1 >= 1, NSC_r2 >= 1) |>
  pull(gene_id)

# Subset GTF to expressed genes, then extract exon only
exons_expressed <- gtf[mcols(gtf)$gene_id %in% expressed_NSC_genes] |>
  subset(type == "exon")

# Collapse overlapping exons per gene to get a non-redundant exon model
exons_by_gene <- split(exons_expressed, exons_expressed$gene_id)
reduced_exons_by_gene <- lapply(exons_by_gene, reduce)

# Compute the gene span per gene
gene_ranges <- lapply(exons_by_gene, function(e) {
  GRanges(
    seqnames = seqnames(e)[1],
    ranges = IRanges(start = min(start(e)), end = max(end(e))),
    strand = strand(e)[1]
  )
})

# Introns = gene span - exons; drop genes with no introns
introns_by_gene <- mapply(setdiff, gene_ranges, reduced_exons_by_gene, SIMPLIFY = FALSE)
introns_by_gene_nonempty <- introns_by_gene[sapply(introns_by_gene, length) > 0]

# Annotate each GRanges entry with the parent gene ID
for (gene in names(introns_by_gene_nonempty)) {
  introns_by_gene_nonempty[[gene]]$gene_id <- gene
}

# Flatten all gene-level GRanges into one data frame
all_introns_df <- names(introns_by_gene_nonempty) |>
  lapply(function(gene_id) {
    gr <- introns_by_gene_nonempty[[gene_id]]
    if (length(gr) == 0) return(NULL)
    data.frame(
      chr = as.character(seqnames(gr)),
      start = start(gr),
      end = end(gr),
      gene_id = gene_id,
      strand = as.character(strand(gr))
    )
  }) |>
  do.call(what = rbind)


# ── 2. Build background intron set ───────────────────────────────────────────

# Load top TIR (transposable-element-associated intronic region) introns
tir_introns <- fread("top_TIR_Neuron_introns_hg38.txt", header = FALSE)
colnames(tir_introns) <- c("chr", "start", "end", "gene_id", "gene_symbol", "strand")

tir_gr <- GRanges(seqnames = tir_introns$chr,
                  ranges   = IRanges(tir_introns$start, tir_introns$end))
all_gr <- GRanges(seqnames = all_introns_df$chr,
                  ranges   = IRanges(all_introns_df$start, all_introns_df$end))

# Remove any all_introns that overlap TIR introns to form the background set
overlaps <- findOverlaps(all_gr, tir_gr)
background_gr <- all_gr[-queryHits(overlaps)]
background_introns <- as.data.frame(background_gr)

# Compute GC content and length for TIR introns
tir_gc <- getSeq(genome, tir_gr) |>
  letterFrequency(letters = c("G", "C"), as.prob = TRUE)
tir_introns <- as_tibble(tir_introns) |>
  mutate(
    GC = rowSums(tir_gc),
    size = end - start + 1,
    index = row_number()
  )

# Compute GC content and length for background introns
bg_gc <- getSeq(genome, background_gr) |>
  letterFrequency(letters = c("G", "C"), as.prob = TRUE)
background_introns <- background_introns |>
  mutate(
    GC = rowSums(bg_gc),
    size = width
  )


# ── 3. Helpful functions ───────────────────────────────────────────────────────

# Match a TIR intron to a background intron with similar size and GC content
match_introns <- function(tir_intron, background_introns,
                          buff_len = 0.1, buff_GC = 0.05) {
  size_range <- tir_intron$size * c(1 - buff_len, 1 + buff_len)
  gc_range <- tir_intron$GC * c(1 - buff_GC,  1 + buff_GC)

  candidates <- background_introns |>
    filter(
      size >= size_range[1], size <= size_range[2],
      GC >= gc_range[1], GC <= gc_range[2]
    )

  if (nrow(candidates) == 0) return(NULL)

  # Pick one matching intron at random
  match <- candidates[sample(nrow(candidates), 1), ]
  match$index <- tir_intron$index
  match
}


# Compare de novo variant (DNV) burden in TIR vs. matched control introns
var_rate_compare <- function(Var_file, Gene_cor_file_target, Gene_cor_file_control) {

  # Overlap variants with an intron set, count per individual ─────────────────
  var_rate <- function(Var_file, Gene_cor_file = NULL) {
    if (!is.null(Gene_cor_file)) {
      Gene_cor <- as.data.table(Gene_cor_file)
      colnames(Gene_cor)[1:3] <- c("CHROM", "POS.start", "POS.end")

      setDT(Var_file)
      setDT(Gene_cor)
      Var_file[, POS.start := POS]
      Var_file[, POS.end := POS]
      setkey(Gene_cor, CHROM, POS.start, POS.end)
      setkey(Var_dt, CHROM, POS.start, POS.end)

      # Overlap-join: retain only variants falling within intron boundaries
      Var_df <- foverlaps(Var_dt, Gene_cor, nomatch = 0L)
    } else {
      Var_df <- Var_file
    }

    # Aggregate DNV count per individual
    dnv_counts <- Var_df |>
      group_by(Ind_ID) |>
      summarise(
        Num_var = n(),
        Sex = first(Sex),
        ASD = first(ASD),
        Fam_ID = first(Fam_ID),
        .groups = "drop"
      )

    # Add zero-count rows for individuals with no overlapping DNVs
    zero_counts <- Var_file |>
      filter(!(Ind_ID %in% dnv_counts$Ind_ID)) |>
      distinct(Ind_ID, Sex, ASD, Fam_ID) |>
      mutate(Num_var = 0L)

    ind_dnv <- bind_rows(dnv_counts, zero_counts) |>
      mutate(
        Num_var = as.numeric(Num_var),
        ASD = as.factor(ASD),
        Sex = as.factor(Sex)
      )

    n_asd <- length(unique(Var_file$Ind_ID[Var_file$ASD == 2]))
    n_nonasd <- length(unique(Var_file$Ind_ID[Var_file$ASD == 1]))

    list(
      asd_rate = sum(Var_df$ASD == 2, na.rm = TRUE) / n_asd,
      nonasd_rate = sum(Var_df$ASD == 1, na.rm = TRUE) / n_nonasd,
      ind_dnv = ind_dnv
    )
  }
  # ───────────────────────────────────────────────────────────────────────────

  target_g <- var_rate(Var_file, Gene_cor_file_target)
  control_g <- var_rate(Var_file, Gene_cor_file_control)

  rates_df <- data.frame(
    row.names = "Rate",
    ASD_target = target_g$asd_rate,
    nonASD_target = target_g$nonasd_rate,
    ASD_control = control_g$asd_rate,
    nonASD_control = control_g$nonasd_rate
  )

  # Combine per-individual counts; sort by family for GEE clustering
  ind_combined <- bind_rows(
    mutate(target_g$ind_dnv, Group = "TIR_introns"),
    mutate(control_g$ind_dnv, Group = "Matched_introns")
  ) |>
    arrange(Fam_ID)

  # Convert p-values to star-notation significance labels
  label_sig <- function(p) {
    case_when(
      p < 0.001 ~ "***",
      p >= 0.001 & p < 0.01 ~ "**",
      p >= 0.01  & p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  }

  # GEE model for ASD probands — identity link to estimate mean difference
  ind_asd <- ind_combined |> filter(ASD == 2)
  gee_asd <- gee(Num_var ~ Group, data = ind_asd, link = "identity", 
                 cond = FALSE, clusterid = ind_asd["Fam_ID"])
  gee_asd_df <- as.data.frame(summary(gee_asd)[["coefficients"]]) |>
    mutate(sig = label_sig(`Pr(>|z|)`), Group = "ASD", variate = c("Matched", "TIR"))

  # GEE model for unaffected siblings
  ind_nonasd <- ind_combined |> filter(ASD == 1)
  gee_nonasd <- gee(Num_var ~ Group, data = ind_nonasd, link = "identity",
                    cond = FALSE, clusterid = ind_nonasd["Fam_ID"])
  gee_nonasd_df <- as.data.frame(summary(gee_nonasd)[["coefficients"]]) |>
    mutate(sig = label_sig(`Pr(>|z|)`), Group = "nonASD", variate = c("Matched", "TIR"))

  list(rates_df, bind_rows(gee_asd_df, gee_nonasd_df))
}


# Run repeatedly sample matched controls and compare DNV rates
run_match_analysis <- function(dnv_data, tir_introns, background_introns, n_iter = 100) {
  rates_list <- vector("list", n_iter)
  gee_list <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    # Sample one matched background intron per TIR intron
    matched_controls <- seq_len(nrow(tir_introns)) |>
      lapply(function(j) match_introns(tir_introns[j, ], background_introns)) |>
      do.call(what = rbind)

    tir_matched <- tir_introns[tir_introns$index %in% matched_controls$index, ]

    result <- var_rate_compare(dnv_data, tir_matched, matched_controls)
    rates_list[[i]] <- result[[1]]
    gee_list[[i]] <- result[[2]]
  }

  rates_combined <- do.call(rbind, rates_list)
  gee_combined <- do.call(rbind, gee_list)

  list(
    rates = as.data.frame(lapply(rates_combined, unlist)),
    gee = as.data.frame(lapply(gee_combined,unlist))
  )
}


# Plot a histogram of 100 times of DNV rates with observed TIR intron rates overlaid
plot_dnv_histogram <- function(rates_df, tir_asd_rate, tir_nonasd_rate,
                               group_colors, cohort_label) {
  plot_df <- data.frame(
    `DNV rate` = c(rates_df$ASD_control, rates_df$nonASD_control),
    Group = rep(c("ASD", "non-ASD"), each = nrow(rates_df)),
    check.names = FALSE
  )

  ggplot(plot_df, aes(x = `DNV rate`, fill = Group)) +
    geom_histogram(position = "identity", alpha = 0.4, bins = 50, color = "black") +
    xlim(0.83, 1.15) +
    # Vertical lines mark the observed TIR intron rates for each group
    geom_vline(aes(xintercept = tir_asd_rate, linetype = "ASD"),
               color = "firebrick3",  linewidth = 1) +
    geom_vline(aes(xintercept = tir_nonasd_rate, linetype = "non-ASD"),
               color = "dodgerblue3", linewidth = 1) +
    scale_fill_manual(values = group_colors) +
    scale_linetype_manual(
      name = "TOP TIR Introns",
      values = c("ASD" = "dashed", "non-ASD" = "dashed")
    ) +
    theme_minimal() +
    labs(
      x = "DNV rate",
      y = "Count",
      fill = "Matched Introns",
      title = paste0(cohort_label,
                     ": Histogram of DNV rate on TOP 55 TIR Introns and matched introns"),
      subtitle = "Background - expressed genes: TPM >= 1 in NSC"
    ) +
    theme(legend.title = element_text(size = 11))
}


# ── 4. Run analyses ───────────────────────────────────────────────

# SPARK cohort — columns: CHROM, POS, REF, ALT, Ind_ID, ASD, Sex, Fam_ID
SPARK_DNV_QC   <- read.delim("SPARK_DNV.tsv")
SPARK_results  <- run_match_analysis(SPARK_DNV_QC, tir_introns, background_introns)
SPARK_gee_TIR  <- SPARK_results$gee |> filter(variate == "TIR")

# SSC cohort — same column structure as SPARK
SSC_DNV_QC   <- read.delim("SSC_DNV.tsv")
SSC_results  <- run_match_analysis(SSC_DNV_QC, tir_introns, background_introns)
SSC_gee_TIR  <- SSC_results$gee |> filter(variate == "TIR")


# ── 5. Visualize DNV rate distributions ──────────────────────────────────────
group_colors <- c(
  "non-ASD" = "dodgerblue3",   
  "ASD" = "firebrick3"
)

SPARK_TOP_matched_intron <- plot_dnv_histogram(
  rates_df = SPARK_results$rates,
  tir_asd_rate = SPARK_results$rates$ASD_target[1],
  tir_nonasd_rate = SPARK_results$rates$nonASD_target[1],
  group_colors = group_colors,
  cohort_label = "SPARK"
)

SSC_TOP_matched_intron <- plot_dnv_histogram(
  rates_df = SSC_results$rates,
  tir_asd_rate = SSC_results$rates$ASD_target[1],
  tir_nonasd_rate = SSC_results$rates$nonASD_target[1],
  group_colors = group_colors,
  cohort_label = "SSC"
)


# ── 6. Significance level distribution ───────────────────────────────────────

# Combine significance counts from both cohorts into one data frame
sig_df <- bind_rows(
  as.data.frame(table(SSC_gee_TIR$Group, SSC_gee_TIR$sig)) |> mutate(Cohort = "SSC"),
  as.data.frame(table(SPARK_gee_TIR$Group, SPARK_gee_TIR$sig)) |> mutate(Cohort = "SPARK")
) |>
  rename(Group = Var1, Sig = Var2, Count = Freq) |>
  mutate(
    Significance = recode(Sig,
      "ns" = "ns (> 0.05)",
      "*" = "0.01-0.05",
      "**" = "0.001-0.01",
      "***" = "< 0.001"
    ),
    Significance = factor(Significance,
                          levels = c("< 0.001", "0.001-0.01", "0.01-0.05", "ns (> 0.05)"))
  )

sig_fig <- 
  ggplot(sig_df, aes(x = Significance, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("ASD" = "firebrick3", "nonASD" = "dodgerblue3")) +
  facet_wrap(~Cohort) +
  theme_bw() +
  labs(
    title = "Significance levels of DNV Rate Differences on TOP 55 TIR Introns vs Matched Introns",
    x = "Significance level of GEE model (p-value)",
    y = "Count",
    fill = "ASD status"
  ) +
  theme(axis.title = element_text(size = 11))
