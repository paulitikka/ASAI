# Reconstruct residue-level RSA from a frozen PISCES chain list and legacy DSSP.
# R >= 4.1; base R only. Run: Rscript reconstruct_pisces_dssp_rsa.R
# Parser check without downloads: Rscript reconstruct_pisces_dssp_rsa.R --self-test
#
# Place the supplied selected_6908_organisms.csv in DATA_DIR. The companion
# aaindex1.txt is used downstream, not for this extraction. For exact paper
# reproduction also supply the 20-AA maximum-ASA table and original DSSP files.
# Downloading current DSSP files produces a new reconstruction; the six analysis
# parts below have the required seven-column layout but will not necessarily
# have the original bundle's hashes, ordering, or counts.

DATA_DIR <- "C:/Users/Pauli Tikka/OneDrive - Karolinska Institutet/Desktop/KI/Lab-specific Work/MasterOrganisms"
PISCES_LIST <- file.path(DATA_DIR, "selected_6908_organisms.csv")
MAX_ASA_CSV <- file.path(DATA_DIR, "max_asa_used_in_paper.csv")
DSSP_DIR <- file.path(DATA_DIR, "dssp_legacy")
DSSP_CHAIN_MAP_CSV <- file.path(DATA_DIR, "dssp_chain_map.csv") # optional; see below
OUTPUT_CSV <- file.path(DATA_DIR, "PISCES_DSSP_RSA_reconstructed.csv")
DOWNLOAD_MISSING_DSSP <- TRUE
WRITE_SIX_ANALYSIS_PARTS <- TRUE # 7-column schema expected by the audit scripts
DSSP_URL <- "https://pdb-redo.eu/dssp/db/%s/legacy"
RSA_CAP <- 1.5
EXPECTED_INPUT_CHAINS <- 6908L
EXPECTED_OUTPUT_CHAINS <- 6342L
EXPECTED_RESIDUES <- 1500829L
AA20 <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]

read_chain_list <- function(path) {
  if (!file.exists(path)) stop("PISCES list missing: ", path)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("pdb", "chain") %in% names(x)))
      names(x)[names(x) == "pdb"] <- "pdb_id"
    if (all(c("pdb_id", "chain") %in% names(x))) {
      out <- x[, c("pdb_id", "chain")]
      out$pdb_id <- tolower(trimws(out$pdb_id))
      out$chain <- trimws(out$chain)
      out$chain[out$chain == ""] <- " "
    } else if ("pdb_chain" %in% names(x)) {
      codes <- x$pdb_chain
      out <- data.frame(pdb_id = tolower(substr(codes, 1L, 4L)),
                        chain = substring(codes, 5L), stringsAsFactors = FALSE)
    } else stop("CSV needs pdb,chain or pdb_id,chain or pdb_chain.")
  } else {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path)
    on.exit(close(con), add = TRUE)
    lines <- trimws(readLines(con, warn = FALSE))
    lines <- lines[nzchar(lines) & !grepl("^#|^PDBchain\\b|^PDB\\s", lines,
                                       ignore.case = TRUE)]
    codes <- vapply(strsplit(lines, "\\s+"), `[`, character(1), 1L)
    if (any(!grepl("^[[:alnum:]]{4}.{1,2}$", codes)))
      stop("Unexpected PISCES first column; supply a pdb_id,chain CSV instead.")
    out <- data.frame(pdb_id = tolower(substr(codes, 1L, 4L)),
                      chain = substring(codes, 5L), stringsAsFactors = FALSE)
  }
  if (anyNA(out) || any(!nzchar(out$pdb_id)) || any(!nzchar(out$chain)) ||
      anyDuplicated(paste(out$pdb_id, out$chain, sep = ":")))
    stop("PISCES chain list has missing or duplicate PDB/chain identifiers.")
  out
}

read_max_asa <- function(path) {
  if (!file.exists(path)) stop("20-AA maximum-ASA table missing: ", path)
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("aa", "max_asa") %in% names(x)))
    stop("Maximum-ASA CSV needs columns aa,max_asa.")
  x$aa <- toupper(trimws(x$aa))
  if (nrow(x) != 20L || anyDuplicated(x$aa) ||
      !setequal(x$aa, AA20) || any(!is.finite(x$max_asa) | x$max_asa <= 0))
    stop("Maximum-ASA CSV must provide one positive finite value per standard AA.")
  setNames(x$max_asa, x$aa)
}

# Legacy DSSP uses one-based columns 6:10 PDB residue number, 11 insertion
# code, 12 chain, 14 AA and 35:38 ACC. The first 5 columns are DSSP's OWN
# serial index and must never replace the PDB residue number in the key.
parse_legacy_dssp <- function(path, pdb_id) {
  lines <- readLines(path, warn = FALSE)
  header <- grep("#\\s+RESIDUE\\s+AA\\s+STRUCTURE", lines)
  if (length(header) != 1L) stop("No unambiguous legacy DSSP residue header.")
  if (header >= length(lines)) stop("No residue records after DSSP header.")
  lines <- lines[seq.int(header + 1L, length(lines))]
  lines <- lines[nchar(lines, type = "bytes") >= 38L]
  if (!length(lines)) stop("No residue records.")
  seqno <- suppressWarnings(as.integer(trimws(substr(lines, 1L, 5L))))
  resseq <- suppressWarnings(as.integer(trimws(substr(lines, 6L, 10L))))
  icode <- substr(lines, 11L, 11L)
  chain <- substr(lines, 12L, 12L)
  aa_raw <- substr(lines, 14L, 14L)
  acc <- suppressWarnings(as.integer(trimws(substr(lines, 35L, 38L))))
  # Lower-case DSSP AA letters denote disulphide-bonded Cys.
  aa <- ifelse(grepl("^[a-z]$", aa_raw), "C", aa_raw)
  keep <- !is.na(seqno) & !is.na(resseq) & aa %in% AA20 &
    !is.na(acc) & acc >= 0L
  data.frame(pdb = toupper(pdb_id), chain = chain[keep],
             dssp_index = seqno[keep], residue_number = resseq[keep],
             insertion_code = icode[keep],
             aa = aa[keep], rsa = rep(NA_real_, sum(keep)),
             ACC = acc[keep], stringsAsFactors = FALSE)
}

read_chain_map <- function(path) {
  if (!file.exists(path))
    return(data.frame(pdb_id = character(), pisces_chain = character(),
                      dssp_chain = character()))
  x <- read.csv(path, stringsAsFactors = FALSE)
  needed <- c("pdb_id", "pisces_chain", "dssp_chain")
  if (!all(needed %in% names(x))) stop("Chain map needs columns: ", paste(needed, collapse = ","))
  x$pdb_id <- tolower(x$pdb_id)
  if (anyDuplicated(paste(x$pdb_id, x$pisces_chain, sep = "|")) ||
      any(nchar(x$dssp_chain) != 1L) || anyNA(x[, needed]))
    stop("Chain map has duplicate, absent or non-legacy DSSP chain IDs.")
  # Two source chains may not point to one legacy DSSP chain in the same PDB.
  if (anyDuplicated(paste(x$pdb_id, x$dssp_chain, sep = "|")))
    stop("Chain map assigns two PISCES chains to one DSSP chain.")
  x[, needed]
}

fetch_dssp <- function(pdb_id, target_dir, download_missing) {
  dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(target_dir, paste0(pdb_id, ".dssp"))
  if (file.exists(path) && file.info(path)$size > 0) return(path)
  if (!download_missing) return(NA_character_)
  url <- sprintf(DSSP_URL, pdb_id)
  tmp <- paste0(path, ".part")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({
    suppressWarnings(utils::download.file(url, tmp, mode = "wb", quiet = TRUE,
                                          method = "libcurl")) == 0L
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tmp) || file.info(tmp)$size < 100L)
    return(NA_character_)
  if (!file.rename(tmp, path)) stop("Could not save downloaded DSSP: ", path)
  path
}

choose_chain <- function(parsed, requested, pdb_id, chain_map) {
  observed <- unique(parsed$chain)
  mapping <- chain_map[chain_map$pdb_id == pdb_id &
                         chain_map$pisces_chain == requested, , drop = FALSE]
  if (nrow(mapping))
    return(if (mapping$dssp_chain %in% observed) mapping$dssp_chain else NA_character_)
  # A multi-character PISCES asym ID cannot be inferred from legacy column 12.
  # Require an explicit, verified author/asym chain mapping in the optional CSV.
  if (nchar(requested) != 1L) return(NA_character_)
  # PISCES sometimes uses 0 for a singleton; check the actual DSSP chains.
  if (requested == "0") {
    if (length(observed) != 1L) return(NA_character_)
    return(observed)
  }
  if (requested %in% observed) requested else NA_character_
}

write_six_parts <- function(full_csv, n_rows) {
  if (n_rows == 0L) return(character())
  con <- file(full_csv, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L)
  per_part <- ceiling(n_rows / 6)
  stem <- sub("\\.csv$", "", full_csv, ignore.case = TRUE)
  paths <- character(6)
  for (i in seq_len(6L)) {
    lines <- readLines(con, n = per_part, warn = FALSE)
    x <- if (length(lines))
      read.csv(text = paste(c(header, lines), collapse = "\n"),
               stringsAsFactors = FALSE, check.names = FALSE) else
      read.csv(text = header, nrows = 0L, check.names = FALSE)
    paths[i] <- paste0(stem, "_part", i, "_of_6.csv")
    write.csv(x[, c("pdb", "chain", "dssp_index", "residue_number",
                    "insertion_code", "aa", "rsa")],
              paths[i], row.names = FALSE, na = "")
  }
  if (length(readLines(con, n = 1L, warn = FALSE)))
    stop("Part splitting omitted rows unexpectedly.")
  paths
}

write_reconstruction <- function(chain_file = PISCES_LIST, max_file = MAX_ASA_CSV,
                                 dssp_dir = DSSP_DIR, output = OUTPUT_CSV,
                                 download_missing = DOWNLOAD_MISSING_DSSP,
                                 chain_map_file = DSSP_CHAIN_MAP_CSV,
                                 compare_expected = TRUE) {
  chains <- read_chain_list(chain_file)
  maximum <- read_max_asa(max_file)
  chain_map <- read_chain_map(chain_map_file)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  tmp_out <- paste0(output, ".part")
  on.exit(unlink(tmp_out), add = TRUE)
  empty <- data.frame(pdb = character(), chain = character(),
                      dssp_index = integer(), residue_number = integer(),
                      insertion_code = character(), aa = character(),
                      rsa = numeric(), ACC = integer())
  write.table(empty, tmp_out, sep = ",", quote = TRUE,
              row.names = FALSE, col.names = TRUE)
  status <- vector("list", nrow(chains))
  grouped <- split(seq_len(nrow(chains)), chains$pdb_id)
  k <- 0L
  for (id in names(grouped)) {
    k <- k + 1L
    if (k %% 100L == 0L) message("PDB entries processed: ", k, "/", length(grouped))
    ix <- grouped[[id]]
    path <- fetch_dssp(id, dssp_dir, download_missing)
    parsed <- if (!is.na(path)) tryCatch(parse_legacy_dssp(path, id),
      error = function(e) structure(list(error = conditionMessage(e)), class = "dssp_error")) else NULL
    for (j in ix) {
      requested <- chains$chain[j]
      chosen <- if (is.data.frame(parsed))
        choose_chain(parsed, requested, id, chain_map) else NA_character_
      msg <- if (is.null(parsed)) "dssp_unavailable" else if
        (inherits(parsed, "dssp_error")) paste0("parse_error: ", parsed$error) else if
        (is.na(chosen) && nchar(requested) > 1L) "multi_char_chain_needs_verified_mapping" else if
        (is.na(chosen)) "chain_not_unambiguous" else "ok"
      part <- if (identical(msg, "ok")) parsed[parsed$chain == chosen, , drop = FALSE] else empty[0, 1:7]
      if (identical(msg, "ok") && nrow(part)) {
        key <- paste(part$pdb, part$chain, part$residue_number,
                     part$insertion_code, sep = "|")
        if (anyDuplicated(key)) stop("Duplicate full residue key in ", id, " chain ", chosen)
        part$insertion_code[part$insertion_code == " "] <- ""
        # For mapped chains retain the original PISCES chain ID in the CSV.
        part$chain <- requested
        part$rsa <- pmin(part$ACC / unname(maximum[part$aa]), RSA_CAP)
        if (any(!is.finite(part$rsa))) stop("Nonfinite RSA in ", id)
        write.table(part, tmp_out, sep = ",", quote = TRUE, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
      } else if (identical(msg, "ok")) msg <- "no_supported_residues"
      status[[j]] <- data.frame(pdb_id = id, requested_chain = requested,
        parsed_chain = if (is.na(chosen)) "" else chosen,
        n_residues = nrow(part), status = msg,
        dssp_file = if (is.na(path)) "" else path,
        dssp_md5 = if (is.na(path)) "" else unname(tools::md5sum(path)),
        stringsAsFactors = FALSE)
    }
  }
  audit <- do.call(rbind, status)
  n_ok <- sum(audit$status == "ok")
  n_res <- sum(audit$n_residues)
  if (!file.copy(tmp_out, output, overwrite = TRUE)) stop("Could not finalize CSV.")
  parts <- if (WRITE_SIX_ANALYSIS_PARTS) write_six_parts(output, n_res) else character()
  audit_file <- paste0(output, ".chain_audit.csv")
  write.csv(audit, audit_file, row.names = FALSE, na = "")
  info <- c(paste("run_time_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE), sep = ""),
    paste0("pisces_list=", normalizePath(chain_file, winslash = "/")),
    paste0("pisces_md5=", unname(tools::md5sum(chain_file))),
    paste0("max_asa_csv=", normalizePath(max_file, winslash = "/")),
    paste0("max_asa_md5=", unname(tools::md5sum(max_file))),
    paste0("chain_map=", if (file.exists(chain_map_file))
      normalizePath(chain_map_file, winslash = "/") else "none"),
    paste0("dssp_download_url=", if (download_missing) DSSP_URL else "none"),
    paste0("rsa_cap=", RSA_CAP), paste0("input_chains=", nrow(chains)),
    paste0("output_chains=", n_ok), paste0("output_residues=", n_res),
    paste0("distinct_PDB_unavailable=", length(unique(audit$pdb_id[audit$status == "dssp_unavailable"]))),
    paste0("analysis_parts=", paste(basename(parts), collapse = ";")),
    paste0("R_version=", R.version.string))
  writeLines(info, paste0(output, ".provenance.txt"))
  if (compare_expected && (nrow(chains) != EXPECTED_INPUT_CHAINS ||
      n_ok != EXPECTED_OUTPUT_CHAINS || n_res != EXPECTED_RESIDUES))
    warning("Counts differ from paper: ", nrow(chains), "/", n_ok, "/", n_res,
            " versus ", EXPECTED_INPUT_CHAINS, "/", EXPECTED_OUTPUT_CHAINS,
            "/", EXPECTED_RESIDUES, ". Inspect chain audit and inputs.")
  message("Written: ", output, " (", n_ok, " chains; ", n_res, " residues)")
  invisible(list(csv = output, parts = parts, audit = audit_file,
                 n_chains = n_ok, n_residues = n_res))
}

self_test <- function() {
  root <- tempfile("dssp_test_"); dir.create(root)
  dir.create(file.path(root, "dssp"))
  list_file <- file.path(root, "list.txt")
  asa_file <- file.path(root, "asa.csv")
  writeLines(c("PDBchain length method", "1abcA 2 XRAY"), list_file)
  write.csv(data.frame(aa = AA20, max_asa = 100), asa_file, row.names = FALSE)
  record <- function(i, res, ins, chain, aa, acc)
    paste0(sprintf("%5d%5d%s%s %s", i, res, ins, chain, aa),
           strrep(" ", 20L), sprintf("%4d", acc))
  writeLines(c("HEADER", "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC",
    record(1L, 10L, "A", "A", "G", 25L),
    record(2L, 10L, "B", "A", "a", 50L)),
    file.path(root, "dssp", "1abc.dssp"))
  result <- write_reconstruction(list_file, asa_file, file.path(root, "dssp"),
    file.path(root, "out.csv"), FALSE, chain_map_file = file.path(root, "no_map.csv"),
    compare_expected = FALSE)
  x <- read.csv(result$csv, stringsAsFactors = FALSE)
  stopifnot(nrow(x) == 2L, identical(x$insertion_code, c("A", "B")),
            identical(x$aa, c("G", "C")),
            isTRUE(all.equal(x$rsa, c(.25, .5))),
            result$n_chains == 1L)
  message("Self-test passed: insertion codes distinguish equal residue numbers.")
}

if (identical(commandArgs(trailingOnly = TRUE), "--self-test")) {
  self_test()
} else if (sys.nframe() == 0L) {
  write_reconstruction()
}
