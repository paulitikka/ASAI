# ASAI stage 1 and follow-up robustness calculations.
# R >= 4.1; base R only. No package installation is needed.
#
# In RStudio: edit the paths below, then Source this entire script.
# Command line: Rscript ASAI_reproduce_stage1.R RSA_INPUT SELECTED_CSV OUTPUT_DIR [AAINDEX_TXT] [AUDIT_ZIP]
# RSA_INPUT: either one reconstructed residue CSV, or a directory containing
# exactly the six residue_RSA_6908_corrected_part*_of_6.csv files.
# Do not use both forms in one run. Original published counts are only a
# comparison target, never a criterion for deleting or selecting records.
# RUN_ADDITIONAL_ONLY = TRUE reads the residue input and produces only new
# chain inventory, cross-fitted PC2, kingdom resampling and chain-level models.
# Set it to FALSE to rerun all original figures/tables plus these new checks.

DATA_DIR <- "C:/Users/Pauli Tikka/OneDrive - Karolinska Institutet/Desktop/KI/Lab-specific Work/MasterOrganisms"
RSA_INPUT <- file.path(DATA_DIR, "residue_RSA_6908_corrected.csv")
# For the newer reconstruction instead, set RSA_INPUT to PISCES_DSSP_RSA_reconstructed.csv.
RUN_ASAI_DISTRIBUTIONS <- TRUE
RUN_ASAI_DISTRIBUTIONS_ONLY <- FALSE # TRUE: produce only the new AA distributions and summaries.
SELECTED_CSV <- file.path(DATA_DIR, "selected_6908_organisms.csv")
OUTPUT_DIR <- file.path(DATA_DIR, "ASAI_stage1_results")
AAINDEX_TXT <- file.path(DATA_DIR, "aaindex1.txt")
AUDIT_ZIP <- file.path(DATA_DIR, "ASAI_audit_bundle_final.zip") # optional for JSON comparisons
REFERENCE_CHAIN_CSV <- "" # optional original corrected RSA CSV or CSV of its pdb,chain IDs
RUN_ADDITIONAL_ONLY <- FALSE
RUN_FIGURES_TABLES_ONLY <- TRUE # Figures 1-3 and their source/supplementary tables only
ADDITIONAL_RESAMPLES <- 120L
PCA_CROSSFIT_SPLITS <- 100L
RUN_EXTENDED <- TRUE  # 20 hold-out splits, 500 permutations, 60+60 PCA stability draws
MODEL_SPLITS <- 20L
PERMUTATIONS <- 500L
STABILITY_DRAWS <- 60L
AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1L]]

write_csv <- function(x, name) write.csv(x, file.path(OUTPUT_DIR, name), row.names = FALSE, na = "")
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

resolve_rsa_files <- function(input) {
  if (!file.exists(input)) stop("RSA input does not exist: ", input)
  if (!dir.exists(input)) return(input)
  p <- list.files(input, "^residue_RSA_6908_corrected_part[1-6]_of_6\\.csv$", full.names = TRUE)
  if (!length(p)) p <- list.files(input, "^PISCES_DSSP_RSA_reconstructed_part[1-6]_of_6\\.csv$", full.names = TRUE)
  if (length(p) != 6L) stop("Expected six matching residue parts in ", input, "; found ", length(p))
  p[order(as.integer(sub(".*_part([1-6])_of_6\\.csv$", "\\1", p)))]
}

read_residues <- function(files) {
  pieces <- lapply(files, function(f) {
    message("Reading ", basename(f))
    x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    needed <- c("pdb", "chain", "residue_number", "insertion_code", "aa", "rsa")
    if (!all(needed %in% names(x))) stop("Missing columns in ", f, ": ",
                                           paste(setdiff(needed, names(x)), collapse = ", "))
    x <- x[, needed]
    x$pdb <- tolower(trimws(as.character(x$pdb)))
    x$chain <- trimws(as.character(x$chain))
    x$insertion_code[is.na(x$insertion_code)] <- ""
    x$insertion_code <- trimws(as.character(x$insertion_code))
    x$aa <- toupper(trimws(as.character(x$aa)))
    x$rsa <- suppressWarnings(as.numeric(x$rsa))
    x$residue_number <- as.character(x$residue_number)
    x
  })
  x <- do.call(rbind, pieces)
  rm(pieces); gc(FALSE)
  if (!nrow(x)) stop("The RSA input has no rows.")
  if (anyNA(x[, c("pdb", "chain", "aa", "rsa")]) ||
      any(!x$aa %in% AA) || any(!is.finite(x$rsa)) ||
      any(x$rsa < 0 | x$rsa > 1.5 + 1e-10))
    stop("Unexpected missing fields, nonstandard amino acids or RSA outside [0, 1.5].")
  if (anyNA(x$residue_number) || any(!nzchar(x$residue_number)))
    stop("Missing PDB residue numbers; insertion-aware duplicate check cannot run.")
  residue_key <- paste(x$pdb, x$chain, x$residue_number, x$insertion_code, sep = "|")
  dup <- duplicated(residue_key)
  bad <- if (any(dup)) x[which(dup)[seq_len(min(sum(dup), 100L))], ] else x[0, ]
  if (nrow(bad)) write_csv(bad, "duplicate_residue_keys_examples.csv")
  n_dup <- sum(dup)
  rm(residue_key, dup, bad); gc(FALSE)
  if (n_dup) warning(n_dup, " repeated PDB+chain+number+insertion keys; inspect report.")
  list(x = x, n_duplicate_keys = n_dup)
}

pool_residues <- function(x) {
  # Retain first-appearance chain order from the residue input.
  chain_key <- paste(x$pdb, x$chain, sep = "|")
  ids <- unique(chain_key)
  first <- match(ids, chain_key)
  aa_ix <- match(x$aa, AA)
  chain_ix <- match(chain_key, ids)
  bin <- (chain_ix - 1L) * 20L + aa_ix
  counts <- tabulate(bin, nbins = length(ids) * 20L)
  # rowsum excludes empty chain/type cells; use its row labels for sparse fill.
  sparse <- rowsum(matrix(x$rsa, ncol = 1L), bin, reorder = FALSE)
  sums <- numeric(length(ids) * 20L)
  sums[as.integer(rownames(sparse))] <- as.vector(sparse[, 1L])
  C <- matrix(counts, nrow = length(ids), byrow = TRUE, dimnames = list(ids, AA))
  S <- matrix(sums, nrow = length(ids), byrow = TRUE, dimnames = list(ids, AA))
  M <- S / C
  M[C == 0L] <- NA_real_
  ztab <- rowsum(matrix(as.integer(x$rsa == 0), ncol = 1L), chain_ix, reorder = FALSE)
  chain_zero <- as.vector(ztab)
  zero_index <- as.integer(rownames(ztab))
  zero <- numeric(length(ids)); zero[zero_index] <- chain_zero
  list(id = ids, pdb = x$pdb[first], chain = x$chain[first],
       L = rowSums(C), cnt = C, sum = S, means = M,
       n_residues = nrow(x), n_zero = sum(x$rsa == 0), zero_per_chain = zero,
       residue_mean = tapply(x$rsa, factor(x$aa, levels = AA), mean))
}

quant_rows <- function(M, idx, Q = 99L) {
  if (!length(idx)) stop("Cannot calculate quantiles for an empty chain group.")
  ans <- t(vapply(seq_along(AA), function(j) {
    v <- M[idx, j]; v <- v[is.finite(v)]
    if (!length(v)) stop("No chain/type means for ", AA[j])
    as.numeric(quantile(v, probs = seq_len(Q)/(Q + 1), type = 7,
                        names = FALSE))
  }, numeric(Q)))
  rownames(ans) <- AA
  ans
}

decompose <- function(rows, mean_rsa, kingdom = FALSE) {
  Q <- ncol(rows)
  centred <- sweep(rows, 2L, colMeans(rows), "-")
  C <- tcrossprod(centred) / Q  # covariance across Q, centred by quantile
  d <- eigen(C, symmetric = TRUE)
  V <- d$vectors
  if (kingdom) {
    av <- (V[seq_along(AA), 1L] + V[length(AA) + seq_along(AA), 1L])
    if (sum(av * mean_rsa) < 0) V[, 1L] <- -V[, 1L]
    k <- match("K", AA); g <- match("G", AA)
    if (mean(V[c(k, k + 20L), 2L]) < mean(V[c(g, g + 20L), 2L])) V[, 2L] <- -V[, 2L]
  } else {
    if (sum(V[, 1L] * mean_rsa) < 0) V[, 1L] <- -V[, 1L]
    if (V[match("K", AA), 2L] < V[match("G", AA), 2L]) V[, 2L] <- -V[, 2L]
    # PC3 has no biological sign convention; align to the printed manuscript.
    if (V[match("M", AA), 3L] < 0) V[, 3L] <- -V[, 3L]
  }
  list(V = V, share = 100 * pmax(d$values, 0) / sum(pmax(d$values, 0)))
}

# Nature-style figure settings. The blue/orange pair is colour-blind-safe and
# is also distinguished by filled versus open symbols. PDF and EPS remain
# vector artwork; TIFF is written at 600 dpi and PNG is a convenient preview.
NATURE_BLUE <- "#0072B2"
NATURE_ORANGE <- "#D55E00"
NATURE_GREY <- "#777777"

# Fixed normalised offsets for Figure 3. The centres are approximately one
# label-circle diameter from their data points in physical plot coordinates.
# Keeping these compact positions in the R source makes the generated artwork
# identical to the reviewed preview instead of allowing device-specific or
# algorithm-specific label placement to differ.
FIGURE3_LABEL_OFFSETS <- data.frame(
  aa = AA,
  dx = c( 0.000000, -0.038044,  0.042000,  0.040338, -0.040338,
          0.042000, -0.036418, -0.021389,  0.034408,  0.005852,
         -0.040338,  0.036418,  0.034408,  0.000000, -0.016700,
         -0.021389, -0.042000,  0.000000, -0.042000, -0.031955),
  dy = c(-0.060000, -0.025420,  0.000000,  0.016709,  0.016709,
          0.000000,  0.029888, -0.051637,  0.034408, -0.059415,
          0.016709, -0.029888,  0.034408,  0.060000,  0.055053,
         -0.051637,  0.000000, -0.060000,  0.000000, -0.038937),
  stringsAsFactors = FALSE
)

nature_par <- function(mar) {
  # The base PostScript device registers Helvetica, not the generic "sans"
  # family. Raster/Cairo devices can retain their portable sans mapping.
  figure_family <- if (identical(names(dev.cur()), "postscript"))
    "Helvetica" else "sans"
  par(family = figure_family, ps = 7, las = 1, bty = "l", lend = "round",
      ljoin = "round", mgp = c(2.15, 0.55, 0), tcl = -0.22,
      xaxs = "r", yaxs = "r", mar = mar,
      cex.axis = 0.95, cex.lab = 1.05, col.axis = "#202020",
      col.lab = "#202020", fg = "#202020")
}

nature_axes <- function(x_at = NULL, x_labels = TRUE, y_at = NULL,
                        y_labels = TRUE) {
  axis(1, at = x_at, labels = x_labels, lwd = 1.35, lwd.ticks = 1.15)
  axis(2, at = y_at, labels = y_labels, lwd = 1.35, lwd.ticks = 1.15)
}

save_nature_figure <- function(stem, draw, width_mm = 183, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  devices <- list(
    pdf = function(path) {
      if (capabilities("cairo"))
        grDevices::cairo_pdf(path, width = w, height = h, family = "sans",
                             pointsize = 7, onefile = TRUE)
      else
        grDevices::pdf(path, width = w, height = h, family = "Helvetica",
                       pointsize = 7, useDingbats = FALSE)
    },
    eps = function(path)
      grDevices::postscript(path, width = w, height = h, horizontal = FALSE,
                            onefile = FALSE, paper = "special",
                            family = "Helvetica", pointsize = 7),
    tiff = function(path)
      grDevices::tiff(path, width = width_mm/25.4, height = height_mm/25.4,
                      units = "in", res = 600, compression = "lzw",
                      family = "sans", pointsize = 7),
    png = function(path)
      grDevices::png(path, width = width_mm/25.4, height = height_mm/25.4,
                     units = "in", res = 300, family = "sans", pointsize = 7)
  )
  for (ext in names(devices)) {
    path <- file.path(OUTPUT_DIR, paste0(stem, ".", ext))
    devices[[ext]](path)
    tryCatch(draw(), finally = dev.off())
  }
}

pc2_dispersion <- function(pool, pc2) {
  # Final2 Methods: one mean RSA for each type in every chain containing it;
  # SD, IQR and range are calculated over those per-chain means. The PCA uses
  # the same means. This is not the SD of the raw residue-level RSA values.
  stats <- do.call(rbind, lapply(seq_along(AA), function(j) {
    v <- pool$means[, j]
    v <- v[is.finite(v)]
    if (length(v) < 2L) stop("Fewer than two chains for ", AA[j])
    data.frame(aa = AA[j], n_chains = length(v), PC2 = pc2[j],
               mean_per_chain_RSA = mean(v), SD_per_chain_mean_RSA = sd(v),
               IQR_per_chain_mean_RSA = IQR(v, type = 7),
               range_per_chain_mean_RSA = diff(range(v)))
  }))
  calc <- function(field, method = "pearson")
    unname(cor(stats$PC2, stats[[field]], method = method))
  summary <- data.frame(
    measure = c("SD Pearson", "SD Spearman", "IQR Pearson", "Range Pearson"),
    r = c(calc("SD_per_chain_mean_RSA"), calc("SD_per_chain_mean_RSA", "spearman"),
          calc("IQR_per_chain_mean_RSA"), calc("range_per_chain_mean_RSA")),
    manuscript_Final2 = c(-.987, -.982, -.869, -.515), n_types = 20L)
  summary$difference <- summary$r - summary$manuscript_Final2
  # Two-sided tests are descriptive: PC2 and dispersion are derived from
  # the same twenty per-type distributions, not independent measurements.
  summary$p_descriptive <- c(cor.test(stats$PC2, stats$SD_per_chain_mean_RSA)$p.value,
    suppressWarnings(cor.test(stats$PC2, stats$SD_per_chain_mean_RSA,
                              method = "spearman", exact = FALSE)$p.value),
    cor.test(stats$PC2, stats$IQR_per_chain_mean_RSA)$p.value,
    cor.test(stats$PC2, stats$range_per_chain_mean_RSA)$p.value)
  # Supplementary Table S3 has two machine-readable parts. If writexl is
  # installed, the same parts are also placed in a single Excel workbook.
  write_csv(stats, "Supplementary_Table_S3a_PC2_dispersion_by_type.csv")
  write_csv(summary, "Supplementary_Table_S3b_PC2_dispersion_correlations.csv")
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list("S3a per type" = stats, "S3b correlations" = summary),
      file.path(OUTPUT_DIR, "Supplementary_Table_S3_PC2_dispersion.xlsx"))
  }

  draw3 <- function() {
    nature_par(c(3.8, 4.0, 0.7, 0.7))
    x <- stats$SD_per_chain_mean_RSA; y <- stats$PC2
    xr <- range(x); yr <- range(y)
    xpad <- diff(xr) * 0.10; ypad <- diff(yr) * 0.10
    plot(x, y, type = "n", xlim = xr + c(-xpad, xpad),
         ylim = yr + c(-ypad, ypad),
         xlab = "SD of per-chain mean RSA (unitless)",
         ylab = "PC2 coordinate (unitless)", axes = FALSE)
    nature_axes()
    fit <- lm(PC2 ~ SD_per_chain_mean_RSA, data = stats)
    abline(fit, col = NATURE_ORANGE, lwd = 2.0)

    # Apply the compact amino-acid-specific offsets used by the distributed
    # preview. Label centres are kept at approximately one circle diameter
    # from their original points while every original point remains visible.
    off <- FIGURE3_LABEL_OFFSETS[match(stats$aa, FIGURE3_LABEL_OFFSETS$aa), ]
    if (anyNA(off$dx) || anyNA(off$dy))
      stop("Missing a Figure 3 label offset for: ",
           paste(stats$aa[is.na(off$dx) | is.na(off$dy)], collapse = ", "))

    # Work in normalised plot coordinates. If an offset points toward or
    # across the fitted line, reflect its normal component. The label then
    # moves farther into the same half-plane as its data point, which makes it
    # geometrically impossible for the straight leader segment to cross the
    # fitted correlation line. This rule applies to all amino-acid types.
    slope_norm <- unname(coef(fit)[2]) * diff(xr) / diff(yr)
    normal <- c(-slope_norm, 1)
    normal <- normal / sqrt(sum(normal^2))
    point_side <- sign((y - predict(fit)) / diff(yr))
    for (i in seq_len(nrow(off))) {
      v <- c(off$dx[i], off$dy[i])
      normal_component <- sum(v * normal)
      side <- point_side[i]
      if (side == 0) side <- if (normal_component < 0) -1 else 1
      if (side * normal_component < 0)
        v <- v - 2 * normal_component * normal
      off$dx[i] <- v[1]
      off$dy[i] <- v[2]
    }
    lab <- data.frame(x = x + off$dx * diff(xr),
                      y = y + off$dy * diff(yr))
    # The circles are painted after the straight leader segments. Thus each
    # visible line reaches the label from the data-point side and stops cleanly
    # at the circle edge instead of continuing beneath the letter.
    segments(x, y, lab$x, lab$y, col = "#666666", lwd = 1.0)
    points(x, y, pch = 16, cex = 0.55, col = "#202020")
    points(lab$x, lab$y, pch = 21, cex = 1.75, lwd = 1.1,
           bg = NATURE_BLUE, col = "#202020")
    text(lab$x, lab$y, stats$aa, col = "white", cex = 0.78,
         font = 2, family = par("family"))
    legend("topright", legend = sprintf("Pearson r = %.3f", summary$r[1]),
           bty = "n", text.font = 2, cex = 0.95)
  }
  save_nature_figure("Figure_3", draw3, width_mm = 136, height_mm = 98)
  list(by_type = stats, correlations = summary)
}

classify_organism <- function(org) {
  prok <- strsplit(paste(
    "escherichia bacillus pseudomonas mycobacterium staphylococcus streptococcus thermus thermotoga",
    "bacteroides vibrio pyrococcus enterococcus helicobacter salmonella klebsiella lactococcus",
    "lactobacillus clostridium corynebacterium rhodobacter synechocystis synechococcus anabaena",
    "nostoc thermococcus archaeoglobus methanocaldococcus methanobacterium halobacterium haloarcula",
    "sulfolobus aeropyrum nanoarchaeum mycoplasma ureaplasma chlamydia rickettsia borrelia treponema",
    "leptospira neisseria haemophilus pasteurella yersinia shigella proteus serratia acinetobacter",
    "burkholderia ralstonia campylobacter deinococcus geobacillus enterobacter citrobacter",
    "stenotrophomonas xanthomonas agrobacterium rhizobium sinorhizobium mesorhizobium bradyrhizobium",
    "caulobacter sphingomonas flavobacterium cytophaga aquifex chlorobium chloroflexus",
    "rhodopseudomonas rhodospirillum thiobacillus nitrosomonas desulfovibrio geobacter shewanella",
    "psychrobacter moraxella halomonas chromohalobacter ferroplasma picrophilus thermoplasma",
    "methanosarcina methanopyrus methanothermobacter methanococcus archaeon thermoproteus",
    "pyrobaculum pyrococcus methanothermus natronomonas haloferax vibrionales streptomyces"), " +")[[1L]]
  euk <- strsplit(paste(
    "homo mus rattus saccharomyces arabidopsis drosophila caenorhabditis bos oryctolagus sus danio",
    "gallus xenopus schizosaccharomyces candida plasmodium trypanosoma leishmania aspergillus",
    "neurospora kluyveromyces pichia komagataella oryza zea nicotiana solanum chlamydomonas",
    "tetrahymena dictyostelium entamoeba giardia trichomonas toxoplasma cryptosporidium spodoptera",
    "trichoplusia bombyx apis strongylocentrotus ciona branchiostoma petromyzon takifugu oryzias",
    "macaca pan canis felis ovis equus cavia mesocricetus cricetulus chlorocebus saimiri mustela",
    "loxodonta ornithorhynchus monodelphis taeniopygia anas meleagris columba physcomitrella",
    "selaginella populus vitis glycine medicago lactuca helianthus sorghum brachypodium physarum",
    "paramecium euglena naegleria cerevisiae pombe"), " +")[[1L]]
  species <- c("homo sapiens", "mus musculus", "rattus norvegicus", "saccharomyces cerevisiae",
    "arabidopsis thaliana", "drosophila melanogaster", "caenorhabditis elegans", "bos taurus",
    "oryctolagus cuniculus", "sus scrofa", "danio rerio", "gallus gallus", "xenopus laevis",
    "schizosaccharomyces pombe", "candida albicans", "plasmodium falciparum", "trypanosoma brucei",
    "leishmania major", "aspergillus niger", "neurospora crassa", "kluyveromyces lactis",
    "pichia pastoris", "komagataella phaffii", "oryza sativa", "zea mays", "nicotiana tabacum",
    "solanum lycopersicum", "chlamydomonas reinhardtii", "tetrahymena thermophila",
    "dictyostelium discoideum", "entamoeba histolytica", "giardia lamblia", "trichomonas vaginalis",
    "toxoplasma gondii", "cryptosporidium parvum", "spodoptera frugiperda", "trichoplusia ni",
    "bombyx mori", "apis mellifera", "strongylocentrotus purpuratus", "ciona intestinalis",
    "branchiostoma lanceolatum", "petromyzon marinus", "takifugu rubripes", "oryzias latipes",
    "macaca mulatta", "pan troglodytes", "canis lupus familiaris", "felis catus", "ovis aries",
    "equus caballus", "cavia porcellus", "mesocricetus auratus", "cricetulus griseus",
    "chlorocebus aethiops", "saimiri sciureus", "mustela putorius furo", "loxodonta africana",
    "ornithorhynchus anatinus", "monodelphis domestica", "taeniopygia guttata",
    "anas platyrhynchos", "meleagris gallopavo", "columba livia")
  o <- tolower(trimws(ifelse(is.na(org), "", org)))
  first <- vapply(strsplit(o, "[ ,;()/]+"), function(x) if (length(x)) x[1L] else "", "")
  out <- rep("unknown", length(o))
  syn <- grepl("synthetic construct", o, fixed = TRUE)
  vir <- !syn & grepl("virus|phage", o)
  out[syn] <- "synthetic"; out[vir] <- "virus"
  left <- out == "unknown" & nzchar(o)
  out[left & o %in% species] <- "euk"
  left <- out == "unknown" & nzchar(o)
  out[left & first %in% prok] <- "prok"
  left <- out == "unknown" & nzchar(o)
  out[left & first %in% euk] <- "euk"
  left <- out == "unknown" & nzchar(o)
  out[left & grepl("\\b(escherichia|bacillus|pseudomonas|mycobacterium|streptococcus|staphylococcus|clostridium|salmonella|vibrio|helicobacter|lactobacillus|lactococcus|mycoplasma|chlamydia|rickettsia|synechocystis|thermus|sulfolobus|halobacterium|methanococcus|archaeon|archaea|bacterium|bacterial)\\b", o, perl = TRUE)] <- "prok"
  left <- out == "unknown" & nzchar(o)
  out[left & grepl("\\b(homo|mus|rattus|saccharomyces|arabidopsis|drosophila|caenorhabditis|danio|xenopus|gallus|bos|sus|oryza|zea|human|mouse|rat|yeast|thale cress|fruit fly|zebrafish)\\b", o, perl = TRUE)] <- "euk"
  out
}

join_selected <- function(pool) {
  if (!file.exists(SELECTED_CSV)) stop("Selected chain input missing: ", SELECTED_CSV)
  x <- read.csv(SELECTED_CSV, stringsAsFactors = FALSE)
  required <- c("pdb", "chain", "sequence", "organism")
  if (!all(required %in% names(x))) stop("Selected chain file needs: ", paste(required, collapse = ", "))
  x$id <- paste(tolower(trimws(x$pdb)), trimws(x$chain), sep = "|")
  if (anyDuplicated(x$id)) stop("Duplicate PDB/chain in selected list.")
  x$kingdom <- classify_organism(x$organism)
  ix <- match(pool$id, x$id)
  if (anyNA(ix)) warning(sum(is.na(ix)), " RSA chains are absent from the selected list.")
  write_csv(data.frame(pdb = pool$pdb, chain = pool$chain, n_residues = pool$L,
                       kingdom = x$kingdom[ix], organism = x$organism[ix],
                       sequence_available = !is.na(x$sequence[ix]) & nzchar(x$sequence[ix])),
            "chain_inclusion_audit.csv")
  list(selected = x, ix = ix, kingdom = x$kingdom[ix],
       sequence = x$sequence[ix])
}

make_divergence <- function(pool) {
  M <- pool$means; present <- is.finite(M)
  lev <- c("exactly 0", paste0(seq(0, 80, 10), "-", seq(10, 90, 10)), "90+")
  by_type <- matrix(0L, nrow = 11L, ncol = 20L, dimnames = list(lev, AA))
  for (j in seq_along(AA)) {
    v <- M[present[, j], j]
    b <- ifelse(v == 0, 1L, pmin(11L, 1L + ceiling(v * 10)))
    by_type[, j] <- tabulate(b, nbins = 11L)
  }
  fractions <- sweep(by_type, 2L, colSums(by_type), "/")
  cv <- apply(fractions, 1L, function(v) sd(v)/mean(v))
  pairs <- rowSums(by_type)
  out <- data.frame(interval = lev, chain_type_pairs = pairs,
                    percentage = 100 * pairs / sum(pairs),
                    cv_across_types = cv,
                    cv_across_types_percent = 100 * cv)
  write_csv(out, "Supplementary_Table_S1.csv")
  write_csv(data.frame(interval = rep(lev, each = 20L), aa = rep(AA, 11L),
                       count = as.vector(t(by_type)), fraction_of_type = as.vector(t(fractions))),
            "Figure_1_source_by_type.csv")
  out
}

joint_result <- function(pool, e, p) {
  if (!length(e) || !length(p)) stop("One of the kingdom groups is empty.")
  r <- rbind(quant_rows(pool$means, e), quant_rows(pool$means, p))
  g <- colSums(pool$sum[c(e, p), , drop = FALSE]) /
       colSums(pool$cnt[c(e, p), , drop = FALSE])
  d <- decompose(r, g, kingdom = TRUE)
  V <- d$V
  t4 <- data.frame(aa = AA, PC1_euk = V[1:20, 1L], PC1_prok = V[21:40, 1L],
                   delta_PC1 = V[1:20, 1L] - V[21:40, 1L],
                   PC2_euk = V[1:20, 2L], PC2_prok = V[21:40, 2L],
                   delta_PC2 = V[1:20, 2L] - V[21:40, 2L])
  list(table = t4, share = d$share,
       t1 = mean(t4$delta_PC1)/(sd(t4$delta_PC1)/sqrt(20)),
       t2 = mean(t4$delta_PC2)/(sd(t4$delta_PC2)/sqrt(20)))
}

window_comparison <- function(pool, e, p) {
  designs <- list("All chains" = c(0, Inf), "150 <= L < 350" = c(150, 350),
                  "200 <= L < 300" = c(200, 300))
  out <- lapply(names(designs), function(label) {
    w <- designs[[label]]
    ei <- e[pool$L[e] >= w[1] & pool$L[e] < w[2]]
    pi <- p[pool$L[p] >= w[1] & pool$L[p] < w[2]]
    d <- joint_result(pool, ei, pi)
    data.frame(comparison = label, n_euk = length(ei), n_prok = length(pi),
               mean_length_euk = mean(pool$L[ei]), mean_length_prok = mean(pool$L[pi]),
               PC1_t_euk_minus_prok = d$t1, PC2_t_euk_minus_prok = d$t2,
               PC1_share_pct = d$share[1], PC2_share_pct = d$share[2],
               PC3_share_pct = d$share[3])
  })
  do.call(rbind, out)
}

figures <- function(s1, t4, share1, share2) {
  draw1 <- function() {
    nature_par(c(4.0, 4.0, 0.7, 0.6))
    xpos <- seq_len(nrow(s1))
    labels <- c("0", "0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4",
                "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8",
                "0.8-0.9", ">0.9")
    cv_pct <- if ("cv_across_types_percent" %in% names(s1))
      s1$cv_across_types_percent else 100 * s1$cv_across_types
    plot(xpos, cv_pct, type = "n", axes = FALSE,
         xlab = "Per-chain mean RSA interval (unitless)",
         ylab = "Coefficient of variability (%)",
         ylim = c(0, max(cv_pct) * 1.08))
    lines(xpos, cv_pct, col = NATURE_BLUE, lwd = 2.25)
    points(xpos, cv_pct, pch = 21, cex = 1.28,
           bg = NATURE_ORANGE, col = "#202020", lwd = 1.15)
    nature_axes(x_at = xpos, x_labels = labels)
  }
  draw2 <- function() {
    old <- par(no.readonly = TRUE); on.exit(par(old))
    layout(matrix(1:2, nrow = 1L), widths = c(1, 1))
    aa_order <- rev(seq_along(AA)); yy <- seq_along(AA)
    panels <- list(
      list(euk = t4$PC1_euk[aa_order], prok = t4$PC1_prok[aa_order],
           panel = "a",
           axis = sprintf("PC1 coordinate (unitless)\nVariance explained: %.2f%%", share2[1])),
      list(euk = t4$PC2_euk[aa_order], prok = t4$PC2_prok[aa_order],
           panel = "b",
           axis = sprintf("PC2 coordinate (unitless)\nVariance explained: %.2f%%", share2[2]))
    )
    for (j in 1:2) {
      z <- panels[[j]]; rr <- range(c(z$euk, z$prok)); pad <- diff(rr) * 0.08
      nature_par(if (j == 1L) c(4.3, 4.0, 2.6, 0.5) else c(4.3, 1.0, 2.6, 0.7))
      # Extra vertical padding keeps the top A and bottom Y symbols and their
      # connecting segments fully inside the plotting region.
      plot(NA_real_, NA_real_, type = "n", axes = FALSE,
           xlim = rr + c(-pad, pad), ylim = c(0.25, length(yy) + 0.75),
           xlab = z$axis, ylab = "", yaxs = "i")
      mtext(z$panel, side = 3, adj = 0, line = 0.45,
            font = 2, cex = 1.14)
      abline(v = 0, col = "#BDBDBD", lwd = 1.0, lty = 3)
      segments(z$prok, yy, z$euk, yy, col = NATURE_GREY, lwd = 1.55)
      points(z$prok, yy, pch = 21, bg = "white", col = NATURE_ORANGE,
             cex = 1.15, lwd = 1.5)
      points(z$euk, yy, pch = 21, bg = NATURE_BLUE, col = "#202020",
             cex = 1.15, lwd = 1.0)
      axis(1, lwd = 1.35, lwd.ticks = 1.15)
      if (j == 1L)
        axis(2, at = yy, labels = AA[aa_order], tick = TRUE,
             lwd = 1.35, lwd.ticks = 1.15, las = 1)
      if (j == 1L)
        mtext("Amino-acid type", side = 2, line = 2.45, las = 0, cex = 0.95)
      if (j == 2L) {
        legend("top", inset = c(0, -0.12), xpd = NA, horiz = TRUE,
               legend = c("Eukaryotes", "Prokaryotes"),
               pch = 21, pt.bg = c(NATURE_BLUE, "white"),
               col = c("#202020", NATURE_ORANGE), pt.cex = 1.15,
               bty = "n", cex = 0.95, text.col = "#202020")
      }
    }
  }
  save_nature_figure("Figure_1", draw1, width_mm = 136, height_mm = 92)
  save_nature_figure("Figure_2", draw2, width_mm = 183, height_mm = 126)
}

chain_predictors <- function(pool, selected) {
  ix <- which(pool$L >= 40L & !is.na(selected$sequence) & nzchar(selected$sequence))
  seqs <- selected$sequence[ix]
  C <- pool$cnt[ix, , drop = FALSE]
  L <- pool$L[ix]
  KD <- c(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-.4,H=-3.2,I=4.5,
          K=-3.9,L=3.8,M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-.8,T=-.7,
          V=4.2,W=-.9,Y=-1.3)
  hydro <- t(vapply(seqs, function(s) {
    a <- strsplit(toupper(s), "", fixed = TRUE)[[1L]]
    h <- unname(KD[a]); h[is.na(h)] <- 0
    w <- if (length(h) >= 11L) max((cumsum(h)[11:length(h)] -
          c(0, head(cumsum(h), -11L)))/11) else 0
    starts <- if (length(a) >= 20L) seq.int(1L, length(a) - 19L, by = 20L) else integer()
    low <- if (length(starts)) mean(vapply(starts, function(i)
      length(unique(a[i:(i + 19L)])) <= 12L, logical(1))) else 0
    c(mean = mean(h), window = w, low = low)
  }, numeric(3)))
  comp <- sweep(C + .5, 1L, L + 20 * .5, "/")
  colnames(comp) <- paste0("comp_", AA)
  X <- data.frame(log10_length = log10(L), comp, check.names = FALSE)
  X$charge <- (-C[, "D"] - C[, "E"] + C[, "K"] + C[, "R"] + .1*C[, "H"])/L
  X$asymmetry_signed <- (C[, "K"] + C[, "R"] + C[, "H"] - C[, "D"] - C[, "E"])/L
  X$asymmetry_absolute <- abs(X$asymmetry_signed)
  X$hydro_mean <- hydro[, "mean"]
  X$hydro_max11 <- hydro[, "window"]
  X$low_complexity <- hydro[, "low"]
  list(idx = ix, X = X)
}

ridge_fit <- function(X, y) {
  A <- crossprod(X)
  lam <- 1e-9 * sum(diag(A))/ncol(X)
  diag(A) <- diag(A) + lam
  tryCatch(as.vector(solve(A, crossprod(X, y))), error = function(e) rep(NA_real_, ncol(X)))
}
r2_test <- function(y, yhat) {
  den <- sum((y - mean(y))^2)
  if (den == 0 || anyNA(yhat)) NA_real_ else 1 - sum((y - yhat)^2)/den
}

context_models <- function(pool, pred) {
  M <- pool$means[pred$idx, , drop = FALSE]
  P <- pred$X
  comp <- paste0("comp_", AA)
  sets <- list(
    length = "log10_length",
    composition = c(comp, "log10_length"),
    charge = c("charge", "asymmetry_signed", "log10_length"),
    hydrophobicity = c("hydro_mean", "hydro_max11", "low_complexity", "log10_length"),
    all = c(comp, "charge", "asymmetry_signed", "hydro_mean", "hydro_max11",
            "low_complexity", "log10_length"),
    all_abs_asymmetry_sensitivity = c(comp, "charge", "asymmetry_absolute",
            "hydro_mean", "hydro_max11", "low_complexity", "log10_length"))
  single <- c("log10_length", "charge", "asymmetry_signed", "hydro_mean",
              "hydro_max11", "low_complexity")
  mats <- lapply(sets, function(fields) cbind(intercept = 1, as.matrix(P[, fields, drop = FALSE])))
  singles <- lapply(single, function(field) cbind(intercept = 1, P[[field]]))
  names(singles) <- single
  corrows <- do.call(rbind, lapply(seq_along(AA), function(j) {
    good <- which(is.finite(M[, j]))
    if (length(good) < 200L) return(NULL)
    do.call(rbind, lapply(single, function(field) {
      data.frame(aa = AA[j], parameter = field, pearson = cor(P[[field]][good], M[good, j]),
                 spearman = cor(P[[field]][good], M[good, j], method = "spearman"))
    }))
  }))
  write_csv(corrows, "context_correlations_by_type.csv")
  set.seed(20260930)
  out <- vector("list", MODEL_SPLITS * 20L)
  k <- 0L
  for (sp in seq_len(MODEL_SPLITS)) {
    for (j in seq_along(AA)) {
      good <- which(is.finite(M[, j]))
      if (length(good) < 200L) next
      perm <- sample(good, length(good), replace = FALSE)
      ntrain <- floor(length(perm)/2)
      tr <- perm[seq_len(ntrain)]; te <- perm[-seq_len(ntrain)]
      ytr <- M[tr, j]; yte <- M[te, j]
      rec <- data.frame(split = sp, aa = AA[j], stringsAsFactors = FALSE)
      for (nm in names(mats)) {
        X <- mats[[nm]]
        b <- ridge_fit(X[tr, , drop = FALSE], ytr)
        rec[[nm]] <- r2_test(yte, as.vector(X[te, , drop = FALSE] %*% b))
      }
      for (nm in names(singles)) {
        X <- singles[[nm]]
        b <- ridge_fit(X[tr, , drop = FALSE], ytr)
        rec[[paste0("single_", nm)]] <- r2_test(yte, as.vector(X[te, , drop = FALSE] %*% b))
      }
      own <- cbind(1, pool$cnt[pred$idx, j]/pool$L[pred$idx])
      b <- ridge_fit(own[tr, , drop = FALSE], ytr)
      rec$single_own_fraction <- r2_test(yte, as.vector(own[te, , drop = FALSE] %*% b))
      k <- k + 1L; out[[k]] <- rec
    }
    message("Hold-out split ", sp, "/", MODEL_SPLITS)
  }
  cells <- do.call(rbind, out[seq_len(k)])
  write_csv(cells, "context_heldout_cells.csv")
  key <- c("length", "composition", "charge", "hydrophobicity", "all", "all_abs_asymmetry_sensitivity")
  table3 <- data.frame(model = key, mean_heldout_R2 = vapply(key, function(z) mean_na(cells[[z]]), 0.0))
  table3$unexplained_percent <- 100 * (1 - table3$mean_heldout_R2)
  write_csv(table3, "Table_3.csv")
  extra <- data.frame(parameter = c(single, "own_fraction"),
                      mean_pearson = c(vapply(single, function(z) mean(corrows$pearson[corrows$parameter == z]), 0.0), NA_real_),
                      mean_spearman = c(vapply(single, function(z) mean(corrows$spearman[corrows$parameter == z]), 0.0), NA_real_),
                      mean_heldout_R2 = vapply(c(paste0("single_", single), "single_own_fraction"),
                                              function(z) mean_na(cells[[z]]), 0.0))
  write_csv(extra, "context_single_parameters.csv")
  list(table3 = table3, correlations = corrows, extra = extra, cells = cells)
}

permutation_correlations <- function(pool, pred, corrows) {
  M <- pool$means[pred$idx, , drop = FALSE]
  P <- pred$X
  fields <- c("log10_length", "charge", "asymmetry_signed", "hydro_mean", "hydro_max11", "low_complexity")
  groups <- lapply(seq_along(AA), function(j) which(is.finite(M[, j])))
  obs <- vapply(fields, function(z) mean(corrows$pearson[corrows$parameter == z]), 0.0)
  obs_abs <- vapply(fields, function(z) mean(abs(corrows$pearson[corrows$parameter == z])), 0.0)
  ge_signed <- ge_abs <- numeric(length(fields))
  set.seed(4242)
  for (i in seq_len(PERMUTATIONS)) {
    shuffles <- lapply(groups, function(g) if (length(g) >= 200) sample.int(length(g)) else integer())
    p <- vapply(seq_along(fields), function(f) {
      z <- fields[f]
      vapply(seq_along(groups), function(j) {
        g <- groups[[j]]
        if (length(g) < 200) return(NA_real_)
        cor(P[[z]][g][shuffles[[j]]], M[g, j])
      }, 0.0)
    }, numeric(20L))
    ge_signed <- ge_signed + (abs(colMeans(p, na.rm = TRUE)) >= abs(obs))
    ge_abs <- ge_abs + (colMeans(abs(p), na.rm = TRUE) >= obs_abs)
    if (i %% 100L == 0L) message("Permutations: ", i, "/", PERMUTATIONS)
  }
  out <- data.frame(parameter = fields, observed_mean_r = obs,
                    p_signed_mean_consistent = (ge_signed + 1)/(PERMUTATIONS + 1),
                    observed_mean_absolute_r = obs_abs,
                    p_mean_absolute_consistent = (ge_abs + 1)/(PERMUTATIONS + 1))
  write_csv(out, "context_permutation_tests.csv")
  out
}

stability <- function(pool, full_v) {
  get_v <- function(idx) decompose(quant_rows(pool$means, idx),
                                  colMeans(pool$means[idx, , drop = FALSE], na.rm = TRUE))$V[, 1:3, drop = FALSE]
  one <- function(type) {
    ans <- matrix(NA_real_, nrow = STABILITY_DRAWS, ncol = 3L)
    for (i in seq_len(STABILITY_DRAWS)) {
      idx <- if (type == "half") sample.int(length(pool$id), floor(length(pool$id)/2))
             else sample.int(length(pool$id), length(pool$id), replace = TRUE)
      v <- get_v(idx)
      ans[i, ] <- vapply(1:3, function(j) abs(cor(v[, j], full_v[, j])), 0.0)
    }
    data.frame(design = type, PC = paste0("PC", 1:3),
               mean_correlation = colMeans(ans), min_correlation = apply(ans, 2, min),
               max_correlation = apply(ans, 2, max))
  }
  set.seed(53)
  rbind(one("half"), one("bootstrap"))
}

aaindex_analysis <- function(path, table1) {
  if (!file.exists(path)) return(NULL)
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  blocks <- strsplit(text, "//", fixed = TRUE)[[1L]]
  old_order <- strsplit("ARNDCQEGHILKMFPSTWYV", "", fixed = TRUE)[[1L]]
  terms <- c("hydrophobicity", "hydropathy", "hydrophilicity", "partition energy",
             "transfer free energy", "solvent accessibility")
  rec <- lapply(blocks, function(b) {
    ln <- strsplit(b, "\n", fixed = TRUE)[[1L]]
    hi <- grep("^\\s*H\\s", ln); ii <- grep("^\\s*I\\s", ln)
    if (!length(hi) || !length(ii)) return(NULL)
    desc <- grep("^\\s*D\\s", ln)
    desc <- desc[desc < ii[1L]]
    d <- paste(trimws(sub("^\\s*D\\s+", "", ln[desc])), collapse = " ")
    nums <- numeric()
    for (k in seq.int(ii[1L], length(ln))) {
      if (k > ii[1L] && grepl("^\\s*[A-Z]\\s", ln[k])) break
      pieces <- regmatches(ln[k], gregexpr("-?[0-9]+(?:\\.[0-9]+)?", ln[k], perl = TRUE))[[1L]]
      nums <- c(nums, as.numeric(pieces))
      if (length(nums) >= 20L) break
    }
    if (length(nums) < 20L) return(NULL)
    v <- nums[seq_len(20L)]; names(v) <- old_order
    list(accession = trimws(sub("^\\s*H\\s+", "", ln[hi[1L]])),
         description = d, hydro = any(vapply(terms, function(t) grepl(t, tolower(d), fixed = TRUE), logical(1))),
         vals = v[AA])
  })
  rec <- Filter(Negate(is.null), rec)
  out <- do.call(rbind, lapply(1:3, function(k) {
    comp <- table1[[paste0("PC", k)]]
    v <- vapply(rec, function(x) abs(cor(x$vals, comp)), 0.0)
    rk <- rank(-v, ties.method = "first")
    data.frame(PC = paste0("PC", k), accession = vapply(rec, `[[`, "", "accession"),
               description = vapply(rec, `[[`, "", "description"),
               hydrophobicity_label = vapply(rec, `[[`, logical(1), "hydro"),
               abs_pearson = v, rank = rk)
  }))
  write_csv(out, "AAindex_correlations.csv")
  aggregate(rank ~ PC, out[out$hydrophobicity_label, ], median)
}

kingdom_controls <- function(pool, e, p) {
  use <- which(pool$L >= 40L)
  f <- pool$cnt[use, , drop = FALSE]/pool$L[use]
  pc <- eigen(crossprod(scale(f, center = TRUE, scale = FALSE))/nrow(f), symmetric = TRUE)$vectors[, 1L]
  score <- as.vector(scale(f, center = TRUE, scale = FALSE) %*% pc)
  all_score <- rep(NA_real_, length(pool$id)); all_score[use] <- score
  e <- e[pool$L[e] >= 40L]; p <- p[pool$L[p] >= 40L]
  pooled <- c(e, p)
  cut_c <- sort(all_score[pooled])[floor(c(1, 2) * length(pooled)/3) + 1L]
  cut_l <- sort(pool$L[pooled])[floor(c(1, 2) * length(pooled)/3) + 1L]
  cb <- function(i) 1L + as.integer(all_score[i] >= cut_c[1]) + as.integer(all_score[i] >= cut_c[2])
  lb <- function(i) 1L + as.integer(pool$L[i] >= cut_l[1]) + as.integer(pool$L[i] >= cut_l[2])
  one <- function(label, ei, pi) {
    if (length(ei) < 60L || length(pi) < 60L) return(data.frame(cell = label,
      n_euk = length(ei), n_prok = length(pi), mean_delta_PC1 = NA_real_, t = NA_real_))
    d <- joint_result(pool, ei, pi)
    data.frame(cell = label, n_euk = length(ei), n_prok = length(pi),
               mean_delta_PC1 = mean(d$table$delta_PC1), t = d$t1)
  }
  comp <- do.call(rbind, lapply(1:3, function(j) one(paste0("C", j), e[cb(e) == j], p[cb(p) == j])))
  len <- do.call(rbind, lapply(1:3, function(j) one(paste0("L", j), e[lb(e) == j], p[lb(p) == j])))
  joint <- do.call(rbind, lapply(1:3, function(j) do.call(rbind,
    lapply(1:3, function(k) one(paste0("C", j, "L", k),
       e[cb(e) == j & lb(e) == k], p[cb(p) == j & lb(p) == k])))))
  x <- rbind(data.frame(control = "composition", comp),
             data.frame(control = "length", len), data.frame(control = "joint", joint))
  write_csv(x, "kingdom_composition_length_controls.csv")
  x
}

length_slopes <- function(pool) {
  out <- do.call(rbind, lapply(seq_along(AA), function(j) {
    idx <- which(pool$L >= 40L & is.finite(pool$means[, j]))
    if (length(idx) < 200L) return(NULL)
    x <- log10(pool$L[idx]); y <- pool$means[idx, j]
    fit <- lm(y ~ x)
    wf <- lm(y ~ x, weights = pool$cnt[idx, j])
    data.frame(aa = AA[j], n_chains = length(idx), slope = coef(fit)[[2L]],
               slope_SE = coef(summary(fit))[2L, 2L], R2 = summary(fit)$r.squared,
               residue_weighted_slope = coef(wf)[[2L]])
  }))
  write_csv(out, "length_slopes_by_type.csv")
  out
}

exact_signflip <- function(d) {
  # All 2^20 combinations, two-sided; zero permutation count has a valid
  # minimum p of 2/2^20 because the observed vector and its negative occur.
  stopifnot(length(d) == 20L, all(is.finite(d)))
  v <- 0
  for (x in d) v <- c(v + x, v - x)
  mean(abs(v) >= abs(sum(d)) - 1e-12)
}

comparison <- function(pool, table1, t3, t4, s1, s2, share1, share2,
                       selected, dispersion) {
  add <- function(metric, observed, manuscript, tol = 0) {
    data.frame(metric = metric, observed = as.numeric(observed), manuscript = as.numeric(manuscript),
               difference = as.numeric(observed) - as.numeric(manuscript), tolerance = tol,
               within_tolerance = abs(as.numeric(observed) - as.numeric(manuscript)) <= tol)
  }
  out <- rbind(add("residue_rows", pool$n_residues, 1500829),
               add("RSA_chains", length(pool$id), 6342),
               add("selected_chains", nrow(selected$selected), 6908),
               add("chain_type_pairs", sum(is.finite(pool$means)), 120822),
               add("kingdom_euk", sum(selected$kingdom == "euk", na.rm = TRUE), 1778),
               add("kingdom_prok", sum(selected$kingdom == "prok", na.rm = TRUE), 2877),
               add("PC1_share_pct", share1[1], 96.95, .005),
               add("PC2_share_pct", share1[2], 2.85, .005),
               add("joint_PC1_share_pct", share2[1], 96.48, .005),
               add("joint_PC2_share_pct", share2[2], 3.23, .005),
               add("PC2_vs_SD_Pearson", dispersion$correlations$r[1], -.987, .00051),
               add("PC2_vs_SD_Spearman", dispersion$correlations$r[2], -.982, .00051),
               add("PC2_vs_IQR_Pearson", dispersion$correlations$r[3], -.869, .00051),
               add("PC2_vs_range_Pearson", dispersion$correlations$r[4], -.515, .00051),
               add("n_chains_C", dispersion$by_type$n_chains[match("C", AA)], 4327),
               add("n_chains_L", dispersion$by_type$n_chains[match("L", AA)], 6321))
  for (j in seq_along(AA)) {
    out <- rbind(out, add(paste0("Table1_meanRSA_", AA[j]), table1$mean_RSA[j],
                          PAPER_T1$mean_RSA[j], .000051),
                      add(paste0("Table1_PC1_", AA[j]), table1$PC1[j], PAPER_T1$PC1[j], .000051),
                      add(paste0("Table4_dPC1_", AA[j]), t4$delta_PC1[j], PAPER_T4$delta_PC1[j], .000101))
  }
  for (i in seq_len(nrow(s1))) {
    out <- rbind(out, add(paste0("S1_pairs_", s1$interval[i]), s1$chain_type_pairs[i],
                          PAPER_S1$chain_type_pairs[i]),
                      add(paste0("S1_CV_", s1$interval[i]), s1$cv_across_types[i],
                          PAPER_S1$cv_across_types[i], .00051))
  }
  for (i in seq_len(nrow(s2))) {
    out <- rbind(out, add(paste0("S2_N_euk_", i), s2$n_euk[i], PAPER_S2$n_euk[i]),
                      add(paste0("S2_N_prok_", i), s2$n_prok[i], PAPER_S2$n_prok[i]))
  }
  if (!is.null(t3)) for (i in seq_len(5L))
    out <- rbind(out, add(paste0("Table3_R2_", t3$model[i]), t3$mean_heldout_R2[i],
                          PAPER_T3$mean_heldout_R2[i], .000051))
  write_csv(out, "manuscript_numeric_comparison.csv")
  out
}

# Reference tables copied as numbers from the supplied ASAI NatComm_Final2.docx.
# Audit JSON can differ from these printed values; never replace computed values with them.
PAPER_T1 <- read.table(text = "aa mean_RSA PC1 PC2 PC3
A 0.2008 -0.0825 0.0746 -0.3895
C 0.0753 -0.3628 0.0758 0.2121
D 0.3557 0.2533 0.1655 0.2498
E 0.3890 0.3258 0.2580 0.2641
F 0.1180 -0.2669 0.1068 0.0368
G 0.2677 0.0593 -0.0400 0.2716
H 0.2634 0.0544 -0.4627 -0.2559
I 0.1145 -0.2743 0.2028 -0.0531
K 0.4169 0.3873 0.1976 0.2718
L 0.1287 -0.2426 0.2999 -0.0130
M 0.1694 -0.1562 -0.5794 0.3913
N 0.3239 0.1867 -0.1572 -0.2966
P 0.3285 0.1955 -0.1423 0.0527
Q 0.3388 0.2194 -0.0376 -0.2612
R 0.3396 0.2200 -0.0172 -0.1120
S 0.2684 0.0645 -0.0053 -0.2155
T 0.2547 0.0346 0.0457 -0.1884
V 0.1319 -0.2352 0.2524 -0.0917
W 0.1371 -0.2254 -0.2216 0.1873
Y 0.1686 -0.1550 -0.0160 -0.0605", header = TRUE)
PAPER_T4 <- read.table(text = "aa delta_PC1
A .0272
C .0342
D .0077
E .0160
F .0253
G .0310
H .0386
I .0275
K .0155
L .0333
M .0340
N .0349
P .0541
Q .0319
R .0361
S .0622
T .0455
V .0347
W .0081
Y .0214", header = TRUE)
PAPER_S1 <- data.frame(chain_type_pairs = c(965,21536,29016,26912,24198,13293,3652,854,261,85,50),
                       cv_across_types = c(2.866,1.149,.682,.604,.820,1.148,1.075,.782,.753,.917,1.558))
PAPER_S2 <- data.frame(n_euk = c(1778,779,348), n_prok = c(2877,1494,741))
PAPER_T3 <- data.frame(mean_heldout_R2 = c(.1164,.1371,.1176,.1168,.1358))

known_reference_conflicts <- function() {
  audit_s1 <- c(965,21535,29017,26913,24196,13295,3651,854,261,85,50)
  x <- data.frame(metric = c("Table3_length_R2", "Table4_PC1_share_pct",
                             "S2_all_PC1_t", "S2_all_PC2_t", paste0("S1_pairs_", seq_len(11))),
                  manuscript = c(.1164, 96.47, -10.14, 5.86, PAPER_S1$chain_type_pairs),
                  bundled_audit = c(.1154, 96.48, 10.14, -5.86, audit_s1))
  x$difference_audit_minus_manuscript <- x$bundled_audit - x$manuscript
  x[x$difference_audit_minus_manuscript != 0, ]
}

final2_narrative_checks <- function() {
  high <- PAPER_T1$aa[order(PAPER_T1$PC2, decreasing = TRUE)[1:5]]
  low <- PAPER_T1$aa[order(PAPER_T1$PC2)[1:5]]
  data.frame(check = c("Largest printed PC2 coordinates",
                       "Smallest printed PC2 coordinates",
                       "Meaning of dispersion in Methods"),
             calculated_or_defined = c(paste(high, collapse = ", "),
                                       paste(low, collapse = ", "),
                                       "between-chain SD of per-chain per-type means"),
             Final2_text = c("calls L,V,I,E,K the lowest PC2 coordinates",
                             "calls M,H,W,N,P the high-dispersion group",
                             "Discussion says residue-to-residue dispersion"))
}

compare_with_audit_json <- function(pool, table1, t3, t4, s1, s2, share1, share2) {
  if (!file.exists(AUDIT_ZIP)) {
    message("Audit ZIP not found; manuscript comparison is still available.")
    return(invisible(NULL))
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("Install jsonlite to compare calculated values against audit JSON; other analyses are complete.")
    return(invisible(NULL))
  }
  read_json <- function(name) {
    con <- unz(AUDIT_ZIP, paste0("results/", name, ".json"))
    on.exit(close(con))
    jsonlite::fromJSON(paste(readLines(con, warn = FALSE), collapse = "\n"))
  }
  d <- read_json("decompose")
  a <- read_json("admissible_parameters")
  w <- read_json("s6_windows")
  add <- function(metric, current, audit, tol) {
    data.frame(metric = metric, calculated = current, audit = audit,
               difference = current - audit, tolerance = tol,
               within_tolerance = abs(current - audit) <= tol)
  }
  t1ref <- d$table1[match(AA, d$table1$t), ]
  t4ref <- d$table2[match(AA, d$table2$t), ]
  out <- rbind(add("n_residues", pool$n_residues, 1500829, 0),
               add("n_chains", length(pool$id), 6342, 0),
               add("Table1_PC1_share", share1[1], d$shares1[1], .005),
               add("Table1_PC2_share", share1[2], d$shares1[2], .005),
               add("Table4_PC1_share", share2[1], d$shares2[1], .005),
               add("Table4_PC2_share", share2[2], d$shares2[2], .005),
               add(paste0("Table1_mean_RSA_", AA), table1$mean_RSA, t1ref$meanRSA, .000051),
               add(paste0("Table1_PC1_", AA), table1$PC1, t1ref$PC1, .000051),
               add(paste0("Table1_PC2_", AA), table1$PC2, t1ref$PC2, .000051),
               add(paste0("Table1_PC3_", AA), table1$PC3, t1ref$PC3, .000051),
               add(paste0("Table4_dPC1_", AA), t4$delta_PC1, t4ref$dPC1, .000101),
               add(paste0("Table4_dPC2_", AA), t4$delta_PC2, t4ref$dPC2, .000101),
               add(paste0("S1_pairs_", s1$interval), s1$chain_type_pairs, d$S5$pairs, 0),
               add(paste0("S1_CV_", s1$interval), s1$cv_across_types, d$S5$cv, .00051))
  model_keys <- c("length", "composition", "charge", "hydrophobicity", "all")
  for (nm in model_keys)
    out <- rbind(out, add(paste0("Table3_R2_", nm),
                          t3$mean_heldout_R2[t3$model == nm], a$models[[nm]]$meanHeldOutR2, .000051))
  for (i in seq_len(nrow(s2)))
    out <- rbind(out, add(paste0("S2_n_euk_", i), s2$n_euk[i], w$nEuk[i], 0),
                      add(paste0("S2_n_prok_", i), s2$n_prok[i], w$nProk[i], 0))
  write_csv(out, "calculation_vs_audit_json.csv")
  invisible(out)
}

chain_membership_inventory <- function(pool, selected) {
  s <- selected$selected
  ix <- match(s$id, pool$id)
  observed <- !is.na(ix)
  inventory <- data.frame(
    pdb = tolower(trimws(s$pdb)), chain = trimws(s$chain),
    organism = s$organism, selected_sequence_length = nchar(s$sequence),
    current_RSA_chain = observed,
    current_RSA_residues = ifelse(observed, pool$L[ix], 0L),
    current_exact_zero_residues = ifelse(observed, pool$zero_per_chain[ix], 0L),
    stringsAsFactors = FALSE)
  write_csv(inventory, "additional_selected_chain_inventory.csv")
  write_csv(inventory[!observed, , drop = FALSE],
            "additional_currently_missing_selected_chains.csv")
  lines <- c(
    paste("Selected chains:", nrow(s)),
    paste("Current RSA chains:", length(pool$id)),
    paste("Current selected chains without RSA:", sum(!observed)),
    paste("Reported corrected RSA chains:", 6342L),
    "Equal chain counts alone cannot identify the corrected dataset's chains.")
  if (nzchar(REFERENCE_CHAIN_CSV)) {
    if (!file.exists(REFERENCE_CHAIN_CSV))
      stop("REFERENCE_CHAIN_CSV not found: ", REFERENCE_CHAIN_CSV)
    ref <- read.csv(REFERENCE_CHAIN_CSV, stringsAsFactors = FALSE,
                    check.names = FALSE)
    if (!"pdb" %in% names(ref) && "pdb_id" %in% names(ref))
      ref$pdb <- ref$pdb_id
    if (!"chain" %in% names(ref) && "requested_chain" %in% names(ref))
      ref$chain <- ref$requested_chain
    if (!all(c("pdb", "chain") %in% names(ref)))
      stop("REFERENCE_CHAIN_CSV needs pdb,chain (or pdb_id,requested_chain).")
    if (anyNA(ref[, c("pdb", "chain")]) ||
        any(!nzchar(trimws(as.character(ref$pdb)))) ||
        any(!nzchar(trimws(as.character(ref$chain)))))
      stop("Reference chain list contains missing or blank identifiers.")
    ref_ids <- unique(paste(tolower(trimws(as.character(ref$pdb))),
                            trimws(as.character(ref$chain)), sep = "|"))
    inventory$reference_chain <- s$id %in% ref_ids
    write_csv(inventory, "additional_reference_membership_comparison.csv")
    lines <- c(lines,
      paste("Reference path:", normalizePath(REFERENCE_CHAIN_CSV)),
      paste("Reference unique chains:", length(ref_ids)),
      paste("Reference chains present in current RSA:", sum(ref_ids %in% pool$id)),
      paste("Current-only chains:", sum(observed & !inventory$reference_chain)),
      paste("Reference-only chains:", sum(!observed & inventory$reference_chain)),
      paste("Reference IDs outside the selected list:", sum(!ref_ids %in% s$id)))
    if (length(ref_ids) != 6342L)
      warning("The optional reference contains ", length(ref_ids),
              " unique chains, rather than the reported 6342.")
  } else {
    lines <- c(lines,
      "No original membership list supplied: the original 566 excluded chains cannot be inferred from these counts.",
      "Optional: set REFERENCE_CHAIN_CSV to an original corrected RSA file or its pdb,chain inventory.")
  }
  writeLines(lines, file.path(OUTPUT_DIR, "additional_membership_summary.txt"))
  invisible(inventory)
}

kingdom_chain_models <- function(pool, selected, e, p) {
  # Each chain is an observation. Cluster SE by PDB in case one structure
  # contributes several chains. This is an association, not a causal estimate.
  use <- c(e[pool$L[e] >= 40L], p[pool$L[p] >= 40L])
  if (length(use) < 100L) stop("Too few classified chains for follow-up models.")
  euk <- as.numeric(selected$kingdom[use] == "euk")
  log_length <- log10(pool$L[use])
  fractions <- pool$cnt[use, , drop = FALSE] / pool$L[use]
  # Drop one fraction so intercept + 19 fractions are identifiable.
  composition <- fractions[, AA[-1L], drop = FALSE]
  Xs <- list(unadjusted = cbind(intercept = 1, euk = euk),
             length = cbind(intercept = 1, euk = euk, log10_length = log_length),
             length_composition = cbind(intercept = 1, euk = euk,
                                         log10_length = log_length, composition))
  fit_one <- function(X, y, pdb) {
    if (nrow(X) <= ncol(X) + 2L || qr(X)$rank != ncol(X))
      return(c(beta_euk = NA_real_, SE_cluster_PDB = NA_real_, p_normal = NA_real_))
    fit <- lm.fit(X, y)
    bread <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
    if (is.null(bread))
      return(c(beta_euk = NA_real_, SE_cluster_PDB = NA_real_, p_normal = NA_real_))
    scores <- rowsum(X * as.vector(fit$residuals), pdb, reorder = FALSE)
    G <- nrow(scores); n <- nrow(X); k <- ncol(X)
    if (G < 2L) stop("Cannot estimate PDB-clustered variance from one PDB.")
    V <- (G/(G - 1)) * ((n - 1)/(n - k)) *
      (bread %*% crossprod(scores) %*% bread)
    se <- sqrt(max(0, V[2L, 2L]))
    beta <- unname(fit$coefficients[2L])
    c(beta_euk = beta, SE_cluster_PDB = se,
      p_normal = if (se > 0) 2 * pnorm(-abs(beta/se)) else NA_real_)
  }
  M <- pool$means[use, , drop = FALSE]
  out <- do.call(rbind, lapply(seq_along(AA), function(j) {
    good <- which(is.finite(M[, j]))
    do.call(rbind, lapply(names(Xs), function(nm) {
      v <- fit_one(Xs[[nm]][good, , drop = FALSE], M[good, j], pool$pdb[use[good]])
      data.frame(aa = AA[j], model = nm,
                 n_euk = sum(euk[good] == 1), n_prok = sum(euk[good] == 0),
                 beta_euk = v["beta_euk"], SE_cluster_PDB = v["SE_cluster_PDB"],
                 CI95_low = v["beta_euk"] - 1.96 * v["SE_cluster_PDB"],
                 CI95_high = v["beta_euk"] + 1.96 * v["SE_cluster_PDB"],
                 p_normal = v["p_normal"], stringsAsFactors = FALSE)
    }))
  }))
  out$p_BH_20_types <- NA_real_
  for (nm in unique(out$model)) {
    ix <- which(out$model == nm)
    out$p_BH_20_types[ix] <- p.adjust(out$p_normal[ix], method = "BH")
  }
  write_csv(out, "additional_kingdom_chain_level_models.csv")
  summary <- do.call(rbind, lapply(unique(out$model), function(nm) {
    x <- out[out$model == nm, ]
    data.frame(model = nm, n_types = nrow(x), n_positive_beta = sum(x$beta_euk > 0, na.rm = TRUE),
               mean_beta = mean(x$beta_euk, na.rm = TRUE),
               n_BH_below_0_05 = sum(x$p_BH_20_types < .05, na.rm = TRUE))
  }))
  write_csv(summary, "additional_kingdom_chain_level_model_summary.csv")
  invisible(out)
}

crossfit_pc2_dispersion <- function(pool) {
  if (PCA_CROSSFIT_SPLITS < 2L) stop("PCA_CROSSFIT_SPLITS must be at least 2.")
  M <- pool$means
  n <- nrow(M)
  set.seed(20261002)
  out <- do.call(rbind, lapply(seq_len(PCA_CROSSFIT_SPLITS), function(b) {
    order <- sample.int(n)
    train <- order[seq_len(floor(n/2))]
    test <- order[-seq_len(floor(n/2))]
    d <- decompose(quant_rows(M, train),
                   colMeans(M[train, , drop = FALSE], na.rm = TRUE))
    test_sd <- apply(M[test, , drop = FALSE], 2L, sd, na.rm = TRUE)
    data.frame(split = b, n_train = length(train), n_test = length(test),
               r_train_PC2_vs_test_SD = unname(cor(d$V[, 2L], test_sd)),
               PC2_share_train_pct = d$share[2L])
  }))
  write_csv(out, "additional_PC2_dispersion_crossfit_draws.csv")
  x <- out$r_train_PC2_vs_test_SD
  summary <- data.frame(n_splits = length(x), mean_r = mean(x), median_r = median(x),
                        q025_r = unname(quantile(x, .025)),
                        q975_r = unname(quantile(x, .975)),
                        n_negative = sum(x < 0))
  write_csv(summary, "additional_PC2_dispersion_crossfit_summary.csv")
  invisible(summary)
}

kingdom_chain_resampling <- function(pool, e, p) {
  if (ADDITIONAL_RESAMPLES < 2L) stop("ADDITIONAL_RESAMPLES must be at least 2.")
  windows <- list(all = c(0, Inf), length_150_350 = c(150, 350),
                  length_200_300 = c(200, 300))
  set.seed(20261003)
  out <- list(); pos <- 0L
  for (nm in names(windows)) {
    w <- windows[[nm]]
    ei <- e[pool$L[e] >= w[1L] & pool$L[e] < w[2L]]
    pi <- p[pool$L[p] >= w[1L] & pool$L[p] < w[2L]]
    if (min(length(ei), length(pi)) < 100L) stop("Too few chains in ", nm)
    for (b in seq_len(ADDITIONAL_RESAMPLES)) {
      ee <- sample(ei, length(ei), replace = TRUE)
      pp <- sample(pi, length(pi), replace = TRUE)
      d <- joint_result(pool, ee, pp)
      pos <- pos + 1L
      out[[pos]] <- data.frame(design = paste0(nm, "_bootstrap"), draw = b,
        n_euk = length(ee), n_prok = length(pp),
        mean_delta_PC1 = mean(d$table$delta_PC1),
        mean_delta_PC2 = mean(d$table$delta_PC2),
        n_positive_delta_PC1 = sum(d$table$delta_PC1 > 0),
        PC1_share_pct = d$share[1L], PC2_share_pct = d$share[2L])
    }
    message("Additional bootstrap ", nm, ": ", ADDITIONAL_RESAMPLES, " draws")
  }
  m <- min(length(e), length(p))
  for (b in seq_len(ADDITIONAL_RESAMPLES)) {
    ee <- sample(e, m, replace = FALSE)
    pp <- sample(p, m, replace = FALSE)
    d <- joint_result(pool, ee, pp)
    pos <- pos + 1L
    out[[pos]] <- data.frame(design = "all_equal_n_subsample", draw = b,
      n_euk = m, n_prok = m,
      mean_delta_PC1 = mean(d$table$delta_PC1),
      mean_delta_PC2 = mean(d$table$delta_PC2),
      n_positive_delta_PC1 = sum(d$table$delta_PC1 > 0),
      PC1_share_pct = d$share[1L], PC2_share_pct = d$share[2L])
  }
  draws <- do.call(rbind, out)
  write_csv(draws, "additional_kingdom_PCA_resampling_draws.csv")
  summarize <- function(x) c(mean = mean(x), q025 = unname(quantile(x, .025)),
                              q975 = unname(quantile(x, .975)))
  summary <- do.call(rbind, lapply(unique(draws$design), function(nm) {
    x <- draws[draws$design == nm, ]
    a <- summarize(x$mean_delta_PC1); b <- summarize(x$mean_delta_PC2)
    data.frame(design = nm, n_draws = nrow(x), n_euk = x$n_euk[1L],
      n_prok = x$n_prok[1L], PC1_mean = a["mean"], PC1_q025 = a["q025"],
      PC1_q975 = a["q975"], PC2_mean = b["mean"], PC2_q025 = b["q025"],
      PC2_q975 = b["q975"], fraction_PC1_positive = mean(x$mean_delta_PC1 > 0),
      fraction_PC2_positive = mean(x$mean_delta_PC2 > 0))
  }))
  write_csv(summary, "additional_kingdom_PCA_resampling_summary.csv")
  invisible(summary)
}

run_additional_checks <- function(pool, selected) {
  chain_membership_inventory(pool, selected)
  e <- which(selected$kingdom == "euk")
  p <- which(selected$kingdom == "prok")
  kingdom_chain_models(pool, selected, e, p)
  crossfit_pc2_dispersion(pool)
  kingdom_chain_resampling(pool, e, p)
  writeLines(c(
    "Additional analysis uses the full current reconstruction unless an original membership CSV is explicitly supplied.",
    "The chain-level models compare eukaryotic and prokaryotic labels in the selected structures; SEs are clustered by PDB.",
    "Composition adjustment uses nineteen of twenty observed amino-acid fractions plus log10 observed chain length.",
    "PCA resampling recomputes both axes within each draw; bootstrap quantiles refer to chain sampling variation.",
    "Equal-n subsampling intervals show subsampling variation, not confidence intervals.",
    "PC2 cross-fitting calculates the PCA in one half and the SD of chain-level means in the other half.",
    "None of these operations infers or forces the original 6342-chain membership."
  ), file.path(OUTPUT_DIR, "additional_checks_readme.txt"))
  message("Additional checks written to ", OUTPUT_DIR)
  invisible(TRUE)
}

# Additional distribution exports, independent of PCA and Low/High selection.
# ASAI: one mean RSA per chain and amino-acid type; chains weighted equally.
# Individual residue RSA is exported separately and is not called chain ASAI.
write_asai_distributions <- function(pool, residues, selected) {
  chain_id <- rep(pool$id, each = length(AA))
  chain_table <- data.frame(pdb = rep(pool$pdb, each = length(AA)),
    chain = rep(pool$chain, each = length(AA)), aa = rep(AA, times = length(pool$id)),
    n_residues = as.vector(t(pool$cnt)), mean_ASAI = as.vector(t(pool$means)))
  chain_table <- chain_table[is.finite(chain_table$mean_ASAI), , drop = FALSE]
  organism <- selected$selected$organism[selected$ix]
  human_ids <- pool$id[!is.na(organism) &
    grepl("^homo sapiens$", trimws(organism), ignore.case = TRUE)]
  chain_table$is_human_annotation <- paste(chain_table$pdb, chain_table$chain,
    sep = "|") %in% human_ids
  write_csv(chain_table, "ASAI_chain_AA_means.csv")
  summarise_types <- function(ct, rt, scope) do.call(rbind, lapply(AA, function(a) {
    v <- ct$mean_ASAI[ct$aa == a]; r <- rt$rsa[rt$aa == a]
    q <- if (length(v)) unname(quantile(v, c(.05,.25,.5,.75,.95), type = 7)) else rep(NA_real_,5)
    data.frame(scope = scope, aa = a, n_chains_with_aa = length(v), n_residues = length(r),
      mean_ASAI = if (length(v)) mean(v) else NA_real_,
      median_ASAI = q[3], SD_chain_ASAI = if (length(v)>1L) sd(v) else NA_real_,
      ASAI_q05 = q[1], ASAI_q25 = q[2], ASAI_q75 = q[4], ASAI_q95 = q[5],
      residue_weighted_mean_RSA = if (length(r)) mean(r) else NA_real_,
      residue_median_RSA = if (length(r)) median(r) else NA_real_)
  }))
  rt_human <- residues[paste(residues$pdb,residues$chain,sep="|") %in% human_ids,,drop=FALSE]
  summary <- rbind(summarise_types(chain_table,residues,"all_chains"),
    summarise_types(chain_table[chain_table$is_human_annotation,,drop=FALSE],rt_human,"human_annotation"))
  write_csv(summary,"ASAI_AA_distribution_summary.csv")
  # Reflected KDE on the documented [0,1.5] RSA range. Reflection reduces
  # artificial negative-RSA density. Empirical quantiles above remain unsmoothed.
  reflected_density <- function(v) {
    if (length(v)<2L) return(NULL)
    bw <- max(stats::bw.nrd0(v), .005)
    z <- stats::density(c(v,-v,3-v),bw=bw,from=0,to=1.5,n=1024L)
    z$y <- 3*z$y; z
  }
  chain_density <- lapply(AA,function(a) reflected_density(chain_table$mean_ASAI[chain_table$aa==a]))
  residue_density <- lapply(AA,function(a) reflected_density(residues$rsa[residues$aa==a]))
  density_table <- do.call(rbind,lapply(seq_along(AA),function(j) {
    do.call(rbind,lapply(c("chain_ASAI","residue_RSA"),function(unit) {
      z <- if (unit=="chain_ASAI") chain_density[[j]] else residue_density[[j]]
      if (is.null(z)) return(NULL)
      data.frame(unit=unit,aa=AA[j],RSA=z$x,density=z$y)
    }))
  }))
  write_csv(density_table,"ASAI_AA_density_curves.csv")
  palette <- rep(c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00",
    "#56B4E9","#000000","#7B61A8","#808000","#8C564B"),2L)
  line_types <- rep(c(1L,2L),each=10L)
  draw <- function() {
    par(mfrow=c(1,2),mar=c(4.6,4.5,3,1),oma=c(4,0,2.5,0),
        las=1,cex=1.05,cex.axis=1,cex.lab=1.12,bty="l",lwd=1.4)
    collections <- list(chain_density,residue_density)
    titles <- c("Per-chain amino-acid means (ASAI)","Individual residue RSA")
    for (panel in 1:2) {
      ds <- collections[[panel]]
      ymax <- max(vapply(ds,function(z) if(is.null(z)) 0 else max(z$y),0.0))
      plot(NA,xlim=c(0,1.5),ylim=c(0,ymax*1.05),xaxs="i",yaxs="i",
        xlab="Relative solvent accessibility (unitless)",
        ylab="Probability density (per RSA unit)",main=titles[panel],cex.main=1.2)
      for (j in seq_along(AA)) if (!is.null(ds[[j]]))
        lines(ds[[j]]$x,ds[[j]]$y,col=palette[j],lty=line_types[j],lwd=2)
    }
    mtext(sprintf("Amino-acid accessibility distributions: %s chains; %s residues",
      format(length(pool$id),big.mark=","),format(nrow(residues),big.mark=",")),
      side=3,outer=TRUE,line=.5,font=2,cex=1.3)
    par(fig=c(0,1,0,.11),new=TRUE,mar=c(0,0,0,0))
    plot.new(); legend("center",legend=AA,col=palette,lty=line_types,lwd=2,
      ncol=10,bty="n",cex=1.15,x.intersp=.6,y.intersp=.85)
  }
  save_nature_figure("ASAI_AA_distributions",draw,width_mm=220,height_mm=125)
  writeLines(c(
    "ASAI distributions use one mean RSA per PDB chain and amino-acid type, with equal weight per chain containing that type.",
    "Individual residue RSA distributions have a different observational unit and are shown separately.",
    "RSA is unitless and retains values above 1 up to the extraction cap of 1.5.",
    "Smooth curves use a boundary-reflected Gaussian density estimate. Summary quantiles use raw unsmoothed values.",
    "Human-only rows are selected by a quote-aware Homo sapiens organism annotation, not by the legacy taxid text column.",
    "The current human subset is not yet an exact UniProt protein inventory.",
    "No length discrimination, AF-tail selection or PCA is required for these distributions."),
    file.path(OUTPUT_DIR,"ASAI_distribution_definitions.txt"))
  invisible(summary)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 5L) stop("Usage: Rscript ASAI_reproduce_stage1.R RSA_INPUT SELECTED_CSV OUTPUT_DIR [AAINDEX_TXT] [AUDIT_ZIP]")
  if (length(args) >= 1L) RSA_INPUT <<- args[1]
  if (length(args) >= 2L) SELECTED_CSV <<- args[2]
  if (length(args) >= 3L) OUTPUT_DIR <<- args[3]
  if (length(args) >= 4L) AAINDEX_TXT <<- args[4]
  if (length(args) >= 5L) AUDIT_ZIP <<- args[5]
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  files <- resolve_rsa_files(RSA_INPUT)
  write_csv(known_reference_conflicts(), "manuscript_vs_supplied_audit_conflicts.csv")
  write_csv(final2_narrative_checks(), "Final2_narrative_checks.csv")
  input <- read_residues(files)
  pool <- pool_residues(input$x)
  selected <- join_selected(pool)
  if (RUN_ASAI_DISTRIBUTIONS) write_asai_distributions(pool, input$x, selected)
  input$x <- NULL; gc(FALSE)
  if (RUN_ASAI_DISTRIBUTIONS_ONLY) return(invisible(TRUE))
  if (RUN_ADDITIONAL_ONLY) {
    run_additional_checks(pool, selected)
    return(invisible(TRUE))
  }
  idx_all <- seq_along(pool$id)
  mean_rsa <- colMeans(pool$means, na.rm = TRUE)
  d1 <- decompose(quant_rows(pool$means, idx_all), mean_rsa)
  table1 <- data.frame(aa = AA, mean_RSA = mean_rsa,
                       PC1 = d1$V[, 1L], PC2 = d1$V[, 2L], PC3 = d1$V[, 3L],
                       residue_weighted_mean_RSA = as.vector(pool$residue_mean[AA]),
                       n_chains_with_aa = colSums(pool$cnt > 0),
                       n_residues_aa = colSums(pool$cnt))
  write_csv(table1, "Table_1.csv")
  dispersion <- pc2_dispersion(pool, table1$PC2)
  s1 <- make_divergence(pool)
  e <- which(selected$kingdom == "euk")
  p <- which(selected$kingdom == "prok")
  d2 <- joint_result(pool, e, p)
  write_csv(d2$table, "Table_4.csv")
  s2 <- window_comparison(pool, e, p)
  write_csv(s2, "Supplementary_Table_S2.csv")
  write_csv(data.frame(PC = c("PC1", "PC2"),
                       mean_shift_euk_minus_prok = c(mean(d2$table$delta_PC1), mean(d2$table$delta_PC2)),
                       paired_t = c(d2$t1, d2$t2),
                       two_sided_exact_signflip_p = c(exact_signflip(d2$table$delta_PC1),
                                                      exact_signflip(d2$table$delta_PC2))),
            "kingdom_exact_signflip.csv")
  write_csv(data.frame(kingdom = c("euk", "prok"), n_chains = c(length(e), length(p)),
                       n_residues = c(sum(pool$L[e]), sum(pool$L[p])),
                       zero_residues = c(sum(pool$zero_per_chain[e]), sum(pool$zero_per_chain[p])),
                       percent_zero_RSA = 100*c(sum(pool$zero_per_chain[e])/sum(pool$L[e]),
                                                sum(pool$zero_per_chain[p])/sum(pool$L[p]))),
            "kingdom_zero_RSA.csv")
  t2 <- data.frame(number = 1:5,
                   parameter = c("Chain length", "Amino-acid composition", "Net charge",
                                 "Charge asymmetry", "Sequence-encoded hydrophobicity"),
                   terms = c("1 (log10)", "20 fractions", "1", "1",
                             "3 (mean, 11-residue maximum, low-complexity fraction)"))
  write_csv(t2, "Table_2.csv")
  figures(s1, d2$table, d1$share, d2$share)
  writeLines(c(
    "FIGURE AND SUPPLEMENTARY TABLE TEXT (additions only)",
    "",
    "Figure 1 | Divergence among amino acid types across the accessibility range. Coefficient of variability, expressed as a percentage, of the twenty type-specific fractions in each interval of per-chain mean RSA. The exact-zero category is separate. Values and observation counts are provided in Supplementary Table S1.",
    "",
    sprintf("Figure 2 | Joint principal coordinates of amino acid accessibility profiles in eukaryotic and prokaryotic chains. Eukaryotic and prokaryotic coordinates from the common forty-row decomposition are shown separately for PC1 (%.2f%%) and PC2 (%.2f%%). Horizontal segments join the two organism groups for the same amino acid type; filled blue circles denote eukaryotes and open orange circles denote prokaryotes. Coordinates are provided in Table 4.", d2$share[1], d2$share[2]),
    "",
    sprintf("Figure 3 | PC2 is associated with between-chain dispersion of amino acid accessibility. Each labelled point represents one amino acid type. Dispersion is the standard deviation, across chains containing that type, of its per-chain mean RSA. The orange line is the ordinary least-squares descriptive fit (Pearson r = %.3f, n = 20 types). Grey leader segments remain in the same half-plane of the fitted line as their data points and do not cross the fitted line. PC2 and dispersion are calculated from the same underlying type distributions; the correlation is therefore descriptive. Values are provided in Supplementary Table S3.", dispersion$correlations$r[1]),
    "",
    "Supplementary Table S3 | PC2 coordinates and dispersion of per-chain mean RSA by amino acid type. Part a reports, for each type, the number of contributing chains, PC2 coordinate, mean, standard deviation, interquartile range and range of the per-chain mean RSA values. Part b reports Pearson or Spearman correlations of PC2 with the dispersion measures. The accompanying p values are descriptive because both quantities are derived from the same twenty type distributions.",
    "",
    "Artwork exports: editable vector PDF and EPS; 600 dpi RGB TIFF; 300 dpi PNG preview. Final widths are 136 mm for Figures 1 and 3 and 183 mm for Figure 2."
  ), file.path(OUTPUT_DIR, "Figure_1_to_3_and_Table_S3_captions.txt"))
  if (RUN_FIGURES_TABLES_ONLY) {
    message("Figures 1-3 and supplementary/source tables written to ", OUTPUT_DIR)
    return(invisible(TRUE))
  }
  slopes <- length_slopes(pool)
  write_csv(data.frame(Q = c(20L, 50L, 99L),
                       PC1_share_pct = c(decompose(quant_rows(pool$means, idx_all, 20L), mean_rsa)$share[1L],
                                         decompose(quant_rows(pool$means, idx_all, 50L), mean_rsa)$share[1L], d1$share[1L]),
                       PC2_share_pct = c(decompose(quant_rows(pool$means, idx_all, 20L), mean_rsa)$share[2L],
                                         decompose(quant_rows(pool$means, idx_all, 50L), mean_rsa)$share[2L], d1$share[2L]),
                       PC3_share_pct = c(decompose(quant_rows(pool$means, idx_all, 20L), mean_rsa)$share[3L],
                                         decompose(quant_rows(pool$means, idx_all, 50L), mean_rsa)$share[3L], d1$share[3L])),
            "quantile_grid_sensitivity.csv")
  pred <- chain_predictors(pool, selected)
  t3 <- context_models(pool, pred)
  if (RUN_EXTENDED) {
    permutation_correlations(pool, pred, t3$correlations)
    write_csv(stability(pool, d1$V), "component_stability.csv")
    write_csv(kingdom_controls(pool, e, p), "kingdom_control_details.csv")
  }
  ar <- aaindex_analysis(AAINDEX_TXT, table1)
  if (!is.null(ar)) write_csv(ar, "AAindex_hydrophobicity_rank_summary.csv")
  check <- comparison(pool, table1, t3$table3, d2$table, s1, s2,
                      d1$share, d2$share, selected, dispersion)
  compare_with_audit_json(pool, table1, t3$table3, d2$table, s1, s2, d1$share, d2$share)
  run_additional_checks(pool, selected)
  write_csv(data.frame(kingdom = c("euk", "prok", "virus", "synthetic", "unknown"),
                       selected_chains = as.numeric(table(factor(selected$selected$kingdom,
                                                    levels = c("euk","prok","virus","synthetic","unknown")))),
                       RSA_chains = as.numeric(table(factor(selected$kingdom,
                                               levels = c("euk","prok","virus","synthetic","unknown"))))),
            "kingdom_counts.csv")
  note <- c(paste("Input files:", paste(normalizePath(files), collapse = "; ")),
            paste("Input MD5:", paste(unname(tools::md5sum(files)), collapse = "; ")),
            paste("RSA rows:", pool$n_residues), paste("RSA chains:", length(pool$id)),
            paste("Duplicate PDB/chain/residue/insertion keys:", input$n_duplicate_keys),
            paste("Exact zero RSA rows:", pool$n_zero),
            paste("Chain-type pairs:", sum(is.finite(pool$means))),
            paste("Sequence-matched chains length >=40:", length(pred$idx)),
            paste("Euk and prok chains:", length(e), length(p)),
            paste("Table 1 PC shares (%):", paste(round(head(d1$share, 4), 4), collapse = ", ")),
            paste("PC2 versus SD of per-chain per-type mean RSA: Pearson",
                  round(dispersion$correlations$r[1], 5), "; Spearman",
                  round(dispersion$correlations$r[2], 5), "; n types = 20"),
            paste("PC2 versus IQR and full range: Pearson",
                  round(dispersion$correlations$r[3], 5), "and",
                  round(dispersion$correlations$r[4], 5)),
            paste("Table 4 PC shares (%):", paste(round(head(d2$share, 3), 4), collapse = ", ")),
            paste("Mean unweighted length slope +/- between-type SD:",
                  round(mean(slopes$slope), 4), "+/-", round(sd(slopes$slope), 4)),
            paste("Manuscript comparisons outside tolerance:", sum(!check$within_tolerance)),
            "Model table labels in the manuscript conceal log10 length: the audit code includes length in every model.",
            "The audit script defines charge asymmetry as signed (positive-negative)/L; the Methods text says absolute.",
            "The audit permutation code compares |mean signed r| with the permuted mean |r|; both consistent tests are reported here.",
            "R uses its own seeded sampling; hold-out and bootstrap draws can differ from the JavaScript audit.",
            "Manuscript and audit reference numbers sometimes disagree; see comparison CSV and the supplied audit results.",
            "The supplied audit JSON predates Final2 and contains no PC2-dispersion calculation; its other comparisons remain available.",
            "Final2 says L,V,I,E,K have the lowest PC2 coordinates; its Table 1 lists them as the five highest.",
            "Final2's per-type SD is across chain-level means; Discussion's 'residue-to-residue' phrase describes a different unit.",
            paste("R:", R.version.string), paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)))
  writeLines(note, file.path(OUTPUT_DIR, "run_summary_and_caveats.txt"))
  message("Finished. See ", OUTPUT_DIR, "; comparisons outside tolerance: ", sum(!check$within_tolerance))
  invisible(list(Table1 = table1, Table3 = t3$table3, Table4 = d2$table,
                 S1 = s1, S2 = s2, PC2_dispersion = dispersion, comparisons = check))
}

main()  # Also runs when the file is sourced in RStudio.
