## Competition CLIP Analysis:
## BedGraph Normalization by Averaging Replicates
## Adapted from CoCLIP pipeline
## Last Edit: 2026-01-12
## Use this after running normalized_bedgraph.sh

library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

baseDir <- "/Volumes/1TB_Data/Competition_CLIP/Luna_4460_Pool_output/normalized_bedGraph/"
bgDir <- paste0(baseDir, "normalized_combined_bedgraph")

setwd(bgDir)
bgFiles <- list.files(bgDir, pattern = "\\.bedgraph$")

for (bgFile in bgFiles) {
  cat("Processing:", bgFile, "\n")

  out_bg <- str_replace(bgFile, ".bedgraph", ".normalized.bedGraph")

  # Parse filename: GROUP.combined.filtered.STRAND.bedgraph
  # e.g., hnRNPC_inRBM25WT.combined.filtered.pos.bedgraph
  parts <- str_split(bgFile, "\\.")[[1]]
  group <- parts[1] # e.g., "hnRNPC_inRBM25WT"
  strand <- parts[4] # e.g., "pos" or "rev"

  # Determine RBP and condition from group name
  if (grepl("hnRNPC", group)) {
    rbp <- "hnRNPC"
    if (grepl("WT", group)) {
      condition <- "RBM25WT"
    } else {
      condition <- "RBM25Mut"
    }
  } else {
    rbp <- "RBM25"
    if (grepl("WT", group)) {
      condition <- "WT"
    } else {
      condition <- "Mut"
    }
  }

  # Format strand for track name
  strand_label <- ifelse(strand == "pos", "(+)", "(-)")

  # Read bedgraph (skip header line from unionbedg)
  bg <- read.delim(bgFile, header = FALSE, skip = 1)

  # Average replicate columns (columns 4 onwards)
  if (ncol(bg) > 4) {
    bg <- bg %>%
      mutate(final = rowMeans(select(., 4:ncol(bg)), na.rm = TRUE))
  } else {
    bg <- bg %>%
      mutate(final = V4)
  }

  # Keep only coordinates and averaged value
  bg <- bg[, c("V1", "V2", "V3", "final")]

  # Create track header
  bg_header <- paste0("track type=bedGraph name=", rbp, "_", condition, "_", strand_label)

  # Write output to same directory as input (normalized_bedGraph base folder)
  outFile <- file.path(baseDir, out_bg)
  writeLines(bg_header, outFile)
  write.table(bg, outFile, row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE, sep = "\t")

  cat("  -> Saved:", out_bg, "\n")
}

cat("\nDone! Normalized bedgraphs saved to:", baseDir, "\n")
