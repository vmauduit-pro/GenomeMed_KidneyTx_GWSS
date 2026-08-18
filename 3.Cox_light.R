#!/usr/bin/env Rscript

usage <- function(status = 0L) {
  cat(
    paste0(
      "Fit per-variant Cox models for donor-recipient genomic mismatches.\n\n",
      "Usage:\n",
      "  Rscript 3.Cox_light.R --pheno FILE --recipient-dosage FILE \\\n",
      "    --donor-dosage FILE --output FILE [options]\n\n",
      "Required arguments:\n",
      "  --pheno FILE              Phenotype/covariate table\n",
      "  --recipient-dosage FILE   Recipient dosage matrix\n",
      "  --donor-dosage FILE       Donor dosage matrix\n",
      "  --output FILE             Output summary-statistics TSV\n\n",
      "Analysis options:\n",
      "  --event NAME              Outcome label [FAILURE]\n",
      "  --event-col NAME          Event indicator column [value of --event]\n",
      "  --time-col NAME           Follow-up column [TIME_<event>]\n",
      "  --pair-id-col NAME        Phenotype ID matching dosage columns [ID]\n",
      "  --hla MODE                epi or none [epi]\n",
      "  --epi-col NAME            Epitopic mismatch column [EPI]\n",
      "  --covariates LIST         Comma-separated model terms; overrides defaults\n",
      "  --legacy-schema MODE      auto, yes, or no [auto]\n",
      "  --cores N                 Parallel Cox workers [1]\n",
      "  --overwrite               Replace an existing output\n",
      "  -h, --help                Show this message\n\n",
      "Default epi covariates:\n",
      "  PC1_R-PC4_R, PC1_D-PC4_D, YEAR, RANK, SEX_R, AGE_R,\n",
      "  tt(AGE_D), TYPE_D, SEX_D, EPI\n\n",
      "Default non-epi covariates:\n",
      "  PC1_R-PC4_R, PC1_D-PC4_D, YEAR, RANK, SEX_R, AGE_R,\n",
      "  TYPE_D, SEX_D, AGE_D, INCOMP_ABDR\n"
    ),
    file = if (status == 0L) stdout() else stderr()
  )
  quit(save = "no", status = status)
}

parse_args <- function(argv) {
  options <- list(
    pheno = NULL,
    `recipient-dosage` = NULL,
    `donor-dosage` = NULL,
    output = NULL,
    event = "FAILURE",
    `event-col` = NULL,
    `time-col` = NULL,
    `pair-id-col` = "ID",
    hla = "epi",
    `epi-col` = "EPI",
    covariates = NULL,
    `legacy-schema` = "auto",
    cores = "1",
    overwrite = FALSE
  )

  if (length(argv) == 0L) usage(1L)
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (token %in% c("-h", "--help")) usage(0L)
    if (token == "--overwrite") {
      options$overwrite <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
    key <- substring(token, 3L)
    if (!key %in% names(options) || key == "overwrite") {
      stop("Unknown option: ", token, call. = FALSE)
    }
    if (i == length(argv) || startsWith(argv[[i + 1L]], "--")) {
      stop("Missing value for ", token, call. = FALSE)
    }
    options[[key]] <- argv[[i + 1L]]
    i <- i + 2L
  }
  options
}

read_table <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  result <- read.table(
    path,
    header = TRUE,
    sep = "",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "#NA", ".")
  )
  if (nrow(result) == 0L) stop(label, " contains no data rows", call. = FALSE)
  if (anyDuplicated(names(result))) stop(label, " has duplicate column names", call. = FALSE)
  result
}

legacy_pheno_names <- c(
  "ID_D", "ID_R", "ID", "FAILURE", "BPAR", "ABMR", "CHR", "TCMR",
  "TIME_FAILURE", "TIME_BPAR", "TIME_ABMR", "TIME_CHR", "TIME_TCMR",
  "YEAR", "RANK", "SEX_R", "AGE_R", "INCOMP_ABDR", "TYPE_D", "SEX_D", "AGE_D",
  paste0("PC", seq_len(20L), "_R"),
  paste0("PC", seq_len(20L), "_D"),
  "EPI"
)

default_covariates <- function(hla, epi_col) {
  pcs <- c(paste0("PC", seq_len(4L), "_R"), paste0("PC", seq_len(4L), "_D"))
  if (hla == "epi") {
    c(pcs, "YEAR", "RANK", "SEX_R", "AGE_R", "tt(AGE_D)", "TYPE_D", "SEX_D", epi_col)
  } else {
    c(pcs, "YEAR", "RANK", "SEX_R", "AGE_R", "TYPE_D", "SEX_D", "AGE_D", "INCOMP_ABDR")
  }
}

parse_covariates <- function(value, hla, epi_col) {
  if (is.null(value)) return(default_covariates(hla, epi_col))
  terms <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  terms[nzchar(terms)]
}

dosage_matrix <- function(data, sample_ids, label) {
  raw <- as.matrix(data[, sample_ids, drop = FALSE])
  numeric_values <- suppressWarnings(as.numeric(raw))
  converted <- matrix(
    numeric_values,
    nrow = nrow(raw),
    ncol = ncol(raw),
    dimnames = dimnames(raw)
  )
  bad <- is.na(converted) & !is.na(raw)
  if (any(bad)) stop(label, " contains non-numeric dosage values", call. = FALSE)
  converted <- round(converted)
  if (any(converted < 0 | converted > 2, na.rm = TRUE)) {
    stop(label, " contains rounded dosages outside 0, 1, or 2", call. = FALSE)
  }
  converted
}

args <- tryCatch(
  parse_args(commandArgs(trailingOnly = TRUE)),
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    usage(1L)
  }
)

required_arguments <- c("pheno", "recipient-dosage", "donor-dosage", "output")
missing_arguments <- required_arguments[vapply(args[required_arguments], is.null, logical(1L))]
if (length(missing_arguments) > 0L) {
  stop(
    "Missing required argument(s): ",
    paste(paste0("--", missing_arguments), collapse = ", "),
    call. = FALSE
  )
}
if (!args$hla %in% c("epi", "none")) stop("--hla must be epi or none", call. = FALSE)
if (!args$`legacy-schema` %in% c("auto", "yes", "no")) {
  stop("--legacy-schema must be auto, yes, or no", call. = FALSE)
}
cores <- suppressWarnings(as.integer(args$cores))
if (is.na(cores) || cores < 1L) stop("--cores must be a positive integer", call. = FALSE)
if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The R package 'survival' is required", call. = FALSE)
}
if (file.exists(args$output) && !args$overwrite) {
  stop("Output exists; use --overwrite to replace it: ", args$output, call. = FALSE)
}

event_col <- if (is.null(args$`event-col`)) args$event else args$`event-col`
time_col <- if (is.null(args$`time-col`)) paste0("TIME_", args$event) else args$`time-col`
pair_id_col <- args$`pair-id-col`
epi_col <- args$`epi-col`
covariates <- parse_covariates(args$covariates, args$hla, epi_col)
if (length(covariates) == 0L) stop("At least one covariate is required", call. = FALSE)

message("Reading phenotype and dosage tables...")
pheno <- read_table(args$pheno, "Phenotype file")

apply_legacy_schema <- args$`legacy-schema` == "yes" ||
  (args$`legacy-schema` == "auto" &&
     ncol(pheno) == length(legacy_pheno_names) &&
     !all(c(pair_id_col, event_col, time_col) %in% names(pheno)))
if (apply_legacy_schema) {
  if (ncol(pheno) != length(legacy_pheno_names)) {
    stop(
      "Legacy phenotype schema requires ", length(legacy_pheno_names),
      " columns; found ", ncol(pheno),
      call. = FALSE
    )
  }
  names(pheno) <- legacy_pheno_names
  message("Applied the legacy 62-column phenotype schema.")
}

formula_variables <- all.vars(as.formula(paste("~", paste(covariates, collapse = " + "))))
required_pheno <- unique(c(pair_id_col, event_col, time_col, formula_variables))
missing_pheno <- setdiff(required_pheno, names(pheno))
if (length(missing_pheno) > 0L) {
  stop("Phenotype column(s) not found: ", paste(missing_pheno, collapse = ", "), call. = FALSE)
}

pheno[[pair_id_col]] <- as.character(pheno[[pair_id_col]])
if (anyNA(pheno[[pair_id_col]]) || any(!nzchar(pheno[[pair_id_col]]))) {
  stop("Pair IDs must be non-missing and non-empty", call. = FALSE)
}
if (anyDuplicated(pheno[[pair_id_col]])) stop("Phenotype pair IDs are duplicated", call. = FALSE)
if (args$hla == "epi") {
  if (!epi_col %in% names(pheno)) stop("EPI column not found: ", epi_col, call. = FALSE)
  keep <- !is.na(pheno[[epi_col]])
  message("EPI filter retained ", sum(keep), " of ", nrow(pheno), " phenotype rows.")
  pheno <- pheno[keep, , drop = FALSE]
}

recipient <- read_table(args$`recipient-dosage`, "Recipient dosage file")
donor <- read_table(args$`donor-dosage`, "Donor dosage file")
info_columns <- c("CHROM", "POS", "ID", "REF", "ALT")
missing_recipient_info <- setdiff(info_columns, names(recipient))
missing_donor_info <- setdiff(info_columns, names(donor))
if (length(missing_recipient_info) > 0L) {
  stop("Recipient dosage file lacks: ", paste(missing_recipient_info, collapse = ", "), call. = FALSE)
}
if (length(missing_donor_info) > 0L) {
  stop("Donor dosage file lacks: ", paste(missing_donor_info, collapse = ", "), call. = FALSE)
}
if (nrow(recipient) != nrow(donor) ||
    !identical(recipient[, info_columns, drop = FALSE], donor[, info_columns, drop = FALSE])) {
  stop("Variant order or CHROM/POS/ID/REF/ALT differs between donor and recipient files", call. = FALSE)
}
variant_info <- recipient[, info_columns, drop = FALSE]

pair_ids <- pheno[[pair_id_col]]
sample_ids <- pair_ids[pair_ids %in% names(recipient) & pair_ids %in% names(donor)]
if (length(sample_ids) == 0L) stop("No phenotype pair IDs match both dosage files", call. = FALSE)
if (length(sample_ids) < length(pair_ids)) {
  message("Using ", length(sample_ids), " of ", length(pair_ids), " phenotype pairs with both dosages.")
}
pheno <- pheno[match(sample_ids, pheno[[pair_id_col]]), , drop = FALSE]

geno_recipient <- dosage_matrix(recipient, sample_ids, "Recipient dosage file")
geno_donor <- dosage_matrix(donor, sample_ids, "Donor dosage file")
rm(recipient, donor)

# Binary donor-derived non-self mismatch. A recipient heterozygote already has
# both alleles, while a homozygous recipient can lack the donor's other allele.
geno_mismatch <- pmin(abs((geno_recipient - geno_donor) * (geno_recipient != 1)), 1)
rm(geno_recipient, geno_donor)

model_data <- pheno
model_data$.TIME <- suppressWarnings(as.numeric(model_data[[time_col]]))
model_data$.EVENT <- suppressWarnings(as.numeric(model_data[[event_col]]))
if (anyNA(model_data$.TIME) || any(!is.finite(model_data$.TIME)) || any(model_data$.TIME < 0)) {
  stop("Follow-up times must be finite, non-missing, and >= 0", call. = FALSE)
}
if (any(grepl("tt(", covariates, fixed = TRUE)) && any(model_data$.TIME <= 0)) {
  stop("Follow-up times must be > 0 when a tt() covariate uses log(time)", call. = FALSE)
}
if (anyNA(model_data$.EVENT) || any(!model_data$.EVENT %in% c(0, 1))) {
  stop("Event indicators must be 0 or 1", call. = FALSE)
}

cox_formula <- as.formula(
  paste0(
    "survival::Surv(.TIME, .EVENT) ~ ",
    paste(c(covariates, ".MISMATCH"), collapse = " + ")
  )
)

fit_one_variant <- function(i) {
  result <- c(COEF = NA_real_, SE.COEF = NA_real_, P = NA_real_, FREQ = NA_real_, FAILED = 1)
  mismatch <- geno_mismatch[i, ]
  if (all(is.na(mismatch))) return(result)
  result[["FREQ"]] <- mean(mismatch > 0, na.rm = TRUE)
  if (length(unique(mismatch[!is.na(mismatch)])) < 2L) return(result)

  variant_data <- model_data
  variant_data$.MISMATCH <- mismatch
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        cox_formula,
        data = variant_data,
        method = "breslow",
        tt = function(x, t, ...) x * log(t),
        model = FALSE,
        x = FALSE,
        y = FALSE
      ),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(result)
  coefficient_index <- match(".MISMATCH", names(fit$coefficients))
  if (is.na(coefficient_index)) return(result)

  beta <- unname(fit$coefficients[[coefficient_index]])
  variance <- unname(fit$var[coefficient_index, coefficient_index])
  if (!is.finite(beta) || !is.finite(variance) || variance <= 0) return(result)
  standard_error <- sqrt(variance)
  result[["COEF"]] <- beta
  result[["SE.COEF"]] <- standard_error
  result[["P"]] <- 2 * pnorm(abs(beta / standard_error), lower.tail = FALSE)
  result[["FAILED"]] <- 0
  result
}

variant_indices <- seq_len(nrow(geno_mismatch))
workers <- min(cores, length(variant_indices))
message(
  "Fitting ", length(variant_indices), " variants across ", length(sample_ids),
  " donor-recipient pairs with ", workers, " worker(s)..."
)
if (workers > 1L && .Platform$OS.type != "windows") {
  fit_results <- parallel::mclapply(variant_indices, fit_one_variant, mc.cores = workers)
} else {
  if (workers > 1L) message("Parallel workers are unavailable on Windows; using one core.")
  fit_results <- lapply(variant_indices, fit_one_variant)
}
fit_results <- do.call(rbind, fit_results)

result <- data.frame(
  variant_info,
  COEF = fit_results[, "COEF"],
  SE.COEF = fit_results[, "SE.COEF"],
  P = fit_results[, "P"],
  FREQ = fit_results[, "FREQ"],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

output_dir <- dirname(normalizePath(args$output, mustWork = FALSE))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  result,
  file = args$output,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t",
  na = "NA",
  eol = "\n"
)

failed <- sum(fit_results[, "FAILED"])
message("Output: ", normalizePath(args$output))
message("Rows:   ", nrow(result))
message("Fits:   ", nrow(result) - failed, " succeeded; ", failed, " failed or were invariant")
message("FREQ is the proportion of pairs carrying the binary mismatch, not allele frequency.")
