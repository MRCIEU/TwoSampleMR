# The one place where opengwasr's data frames become TwoSampleMR's.
#
# Everything that calls opengwasr, and everything that maps its columns to TwoSampleMR's formats,
# is here: the public functions (extract_instruments(), extract_outcome_data(), available_outcomes(),
# add_metadata(), mv_extract_exposures()) check their arguments and call these.

# --- Functions that have not moved yet -------------------------------------------------------------

# The old OpenGWAS client is gone, and these functions have no opengwasr equivalent yet. They stop, with an error of
# class `twosamplemr_not_yet_available` that names the opengwasr ticket that will bring them back.
not_yet_available <- function(fn, needs, ticket) {
  stop(structure(
    class = c("twosamplemr_not_yet_available", "error", "condition"),
    list(
      message = paste0(
        "`",
        fn,
        "()` is not available yet: it needs ",
        needs,
        ", which opengwasr does not have. See ",
        ticket,
        "."
      ),
      call = NULL,
      fn = fn,
      ticket = ticket
    )
  ))
}

# --- Messages -------------------------------------------------------------------------------------

# A message of class `twosamplemr_deprecated_argument`, so that a script (or a test) can tell it from
# any other message.
deprecated_argument_message <- function(text) {
  message(structure(
    class = c("twosamplemr_deprecated_argument", "message", "condition"),
    list(message = paste0(text, "\n"), call = NULL)
  ))
}

# `opengwas_jwt` is read from the Source now.
deprecate_opengwas_jwt <- function(fn) {
  deprecated_argument_message(paste0(
    "`opengwas_jwt` is no longer used by `",
    fn,
    "()` and is ignored. ",
    "Credentials come from the Source: `opengwasr::og_api()` reads `OPENGWAS_JWT`, ",
    "or pass your own Source with `src`."
  ))
}

# The proxy arguments of extract_outcome_data() and mv_extract_exposures(), named by whichever ones
# the user passed.
deprecate_proxy_arguments <- function(passed, fn = "extract_outcome_data") {
  deprecated_argument_message(paste0(
    "The ",
    paste0("`", passed, "`", collapse = ", "),
    " argument",
    if (length(passed) > 1) "s" else "",
    " of `",
    fn,
    "()` no longer ",
    if (length(passed) > 1) "have" else "has",
    " any effect. LD proxies are not looked up: the Reference-Completed Store Releases that ",
    "opengwasr reads already hold Imputed Associations for most variants a study did not report. ",
    "They come back as the variant itself, with `association_status.outcome` set to \"imputed\" ",
    "(see `harmonise_data()` for the count)."
  ))
}

deprecate_splitsize <- function() {
  deprecated_argument_message(paste0(
    "`splitsize` is no longer needed by `extract_outcome_data()` and is ignored: ",
    "opengwasr splits large queries itself."
  ))
}

# --- Analyses -------------------------------------------------------------------------------------

# Identifiers as the old OpenGWAS client wrote them: a bare number is an `ieu-a-` id (or `met-a-`,
# `met-b-`, `met-c-` in the ranges those batches took over), and `UKB-a:` style prefixes are
# `ukb-a-`. Anything else is passed to opengwasr as it is.
legacy_analysis_ids <- function(x) {
  x <- as.character(x)
  bare <- !is.na(x) & !grepl("\\D", x)
  x[bare] <- paste0("ieu-a-", x[bare])
  for (prefix in c("UKB-a:", "UKB-b:", "UKB-c:", "IEU-a:")) {
    index <- startsWith(x, prefix)
    new <- paste0(tolower(sub(":", "", prefix)), "-")
    x[index] <- sub(prefix, new, x[index], fixed = TRUE)
  }
  to_met <- function(x, range, batch) {
    index <- x %in% paste0("ieu-a-", range)
    x[index] <- sub("ieu-a-", batch, x[index], fixed = TRUE)
    x
  }
  x <- to_met(x, 303:754, "met-a-")
  x <- to_met(x, 119:269, "met-b-")
  to_met(x, 838:960, "met-c-")
}

# The canonical ids of the Analyses `ids` name. An id that names no Analysis is warned about and
# left out (as TwoSampleMR always did), unless none is known, which is the opengwasr error.
resolve_analyses <- function(src, ids) {
  ids <- unique(legacy_analysis_ids(ids))
  info <- tryCatch(
    opengwasr::og_analyses(src, ids),
    opengwasr_analysis_not_found = function(e) {
      known <- setdiff(ids, e$ids)
      if (length(known) == 0) {
        stop(e)
      }
      warning(structure(
        class = c("twosamplemr_unknown_analyses", "warning", "condition"),
        list(
          message = paste0(
            "Left out because they name no Analysis: ",
            paste(e$ids, collapse = ", "),
            "."
          ),
          call = NULL,
          ids = e$ids
        )
      ))
      opengwasr::og_analyses(src, known)
    }
  )
  info$id
}

# The catalogue, as `available_outcomes()` returns it.
catalogue <- function(src) {
  a <- opengwasr::og_analyses(src)
  a$trait <- a$analysis_label
  a
}

# The catalogue for the functions that merge study characteristics into MR results. The old
# catalogue's `subcategory` and `nsnp` are not in the new one: results are not grouped by
# subcategory, and the number of SNPs on the array is not reported.
catalogue_for_results <- function() {
  ao <- available_outcomes()
  if (!"subcategory" %in% names(ao)) {
    message(
      "The OpenGWAS catalogue no longer has `subcategory` or `nsnp`, so outcomes are not ",
      "grouped by subcategory and the number of SNPs on the array is not reported."
    )
    ao$subcategory <- ""
  }
  ao
}

# --- Exposure data --------------------------------------------------------------------------------

# Top hits of `analyses` (canonical ids) in TwoSampleMR's exposure format.
opengwasr_instruments <- function(src, analyses, p1, clump, r2, kb) {
  d <- opengwasr::og_tophits(
    src,
    analyses,
    threshold = p1,
    limit = 10000,
    clump = as.logical(clump),
    r2 = r2,
    kb = kb
  )
  if (nrow(d) == 0) {
    return(NULL)
  }
  d$SNP <- variant_name(d)
  # Case and control counts of a continuous trait are all NA, and not worth a column.
  for (column in c("ncase", "ncontrol")) {
    if (all(is.na(d[[column]]))) d[[column]] <- NULL
  }
  d$phenotype <- paste0(d$analysis_label, " || id:", d$id)
  key <- paste(d$id, tolower(d$SNP))
  res <- format_data(
    d,
    type = "exposure",
    snps = NULL,
    phenotype_col = "phenotype",
    snp_col = "SNP",
    chr_col = "chromosome",
    pos_col = "position",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "n",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    min_pval = 1e-200,
    id_col = "id"
  )
  index <- match(paste(res$id.exposure, res$SNP), key)
  # format_data() lower-cases SNP, which would break an ALID.
  alid <- is.na(d$rsid[index])
  res$SNP[alid] <- d$alid[index][alid]
  res <- carry_columns(
    res,
    d,
    index,
    c("effect_scale", "association_status", "clump_ld", "clump_ld_choice", "clump_r2", "clump_kb"),
    "exposure"
  )
  res$data_source.exposure <- "igd"
  res
}

# --- Outcome data ---------------------------------------------------------------------------------

# Associations of `snps` in `analyses` (canonical ids) in TwoSampleMR's outcome format.
opengwasr_outcomes <- function(src, snps, analyses) {
  d <- opengwasr::og_associations(src, snps, analyses)
  if (is.null(nrow(d)) || nrow(d) == 0) {
    return(NULL)
  }
  d <- format_d(d, snps)
  if (nrow(d) > 0) {
    d$data_source.outcome <- "igd"
    return(d)
  }
  NULL
}

#' Format the table that opengwasr returns
#'
#' @param d Data frame from [opengwasr::og_associations()]
#' @param snps The SNPs that were asked for, so that a variant asked for by ALID keeps its ALID.
#' @return Data frame
#' @keywords internal
format_d <- function(d, snps = NULL) {
  d1 <- data.frame(
    SNP = variant_name(d, snps),
    chr = as.character(d$chromosome),
    pos = as.character(d$position),
    beta.outcome = as.numeric(d$beta),
    se.outcome = as.numeric(d$se),
    samplesize.outcome = as.numeric(d$n),
    pval.outcome = as.numeric(d$pval),
    eaf.outcome = as.numeric(d$eaf),
    effect_allele.outcome = as.character(d$effect_allele),
    other_allele.outcome = as.character(d$other_allele),
    outcome = as.character(d$analysis_label),
    id.outcome = as.character(d$id),
    stringsAsFactors = FALSE
  )
  if (!all(is.na(d$ncase))) {
    d1$ncase.outcome <- as.numeric(d$ncase)
  }
  if (!all(is.na(d$ncontrol))) {
    d1$ncontrol.outcome <- as.numeric(d$ncontrol)
  }
  d1 <- carry_columns(d1, d, seq_len(nrow(d)), c("effect_scale", "association_status"), "outcome")
  d <- d1

  if (nrow(d) == 0) {
    message("No matches")
    return(d)
  }

  d$originalname.outcome <- d$outcome
  d$outcome.deprecated <- paste0(d$outcome, " || ", d$consortium.outcome, " || ", d$year.outcome)
  d$outcome <- paste0(d$outcome, " || id:", d$id.outcome)

  rem <- is.na(d$beta.outcome) & is.na(d$pval.outcome)
  d <- subset(d, !rem)

  # For any that have missing SE but available beta and pval, infer SE
  index <- (is.na(d$se.outcome) | d$se.outcome == 0) &
    (!is.na(d$beta.outcome) & !is.na(d$pval.outcome))
  if (any(index)) {
    d$se.outcome[index] <- get_se(d$beta.outcome[index], d$pval.outcome[index])
  }

  d <- cleanup_outcome_data(d)

  mrcols <- c("beta.outcome", "se.outcome", "effect_allele.outcome")
  d$mr_keep.outcome <- complete.cases(d[, mrcols])
  if (!all(d$mr_keep.outcome)) {
    missinginfosnps <- paste(subset(d, !mr_keep.outcome)$SNP, collapse = " ")
    warning(
      "The following SNP(s) are missing required information for the MR tests and will be excluded: ",
      missinginfosnps
    )
  }
  if (!any(d$mr_keep.outcome)) {
    warning(
      "None of the provided SNPs can be used for MR analysis, they are missing required information."
    )
  }

  return(d)
}

# --- Metadata -------------------------------------------------------------------------------------

# What each metadata field of add_metadata() is called in `opengwasr::og_analyses()`.
analyses_names <- c(
  sample_size = "sample_size",
  ncase = "n_cases",
  ncontrol = "n_controls",
  effect_scale = "stored_effect_scale",
  sd = "original_sd"
)

# The `units` label that effect scales had before: `units.exposure` and `units.outcome` are read by
# add_rsq(), standardise_units(), get_rsq() and run_mrmix().
units_from_effect_scale <- function(effect_scale) {
  units <- c(sd = "SD", log_or = "log odds")[effect_scale]
  names(units) <- NULL
  unknown <- unique(effect_scale[is.na(units) & !is.na(effect_scale)])
  if (length(unknown) > 0) {
    message(
      "There is no `units` label for the effect scale ",
      paste0("\"", unknown, "\"", collapse = ", "),
      ", so it is NA."
    )
  }
  units
}

# Metadata of the Analyses `id` (as the data name them, which may be legacy ids), with columns named
# `<field>.<what>` as add_metadata() adds them. NULL where none of the ids is an Analysis.
opengwasr_metadata <- function(src, id, what, cols) {
  info <- tryCatch(
    opengwasr::og_analyses(src, id),
    opengwasr_analysis_not_found = function(e) {
      known <- setdiff(id, e$ids)
      if (length(known) == 0) data.frame() else opengwasr::og_analyses(src, known)
    }
  )
  if (nrow(info) == 0) {
    message(what, ": none of the IDs found in database")
    return(NULL)
  }
  # Match on the id the data carry, which may be a legacy id.
  info$id <- info$requested_id
  for (col in cols) {
    from <- if (col %in% names(analyses_names)) analyses_names[[col]] else col
    info[[col]] <- if (from %in% names(info)) info[[from]] else NA
  }
  if ("effect_scale" %in% cols) {
    info$units <- units_from_effect_scale(info$effect_scale)
    cols <- c(cols, "units")
  }
  info <- subset(info, select = c("id", cols))
  names(info) <- paste0(names(info), ".", what)
  names(info)[names(info) == paste0("sample_size.", what)] <- paste0("samplesize.", what)
  info
}

# `||` for NULL, which base R has only from 4.4.
`%||%` <- function(x, y) if (is.null(x)) y else x

# A Local Source (`opengwasr::og_local()`): Store Releases on this machine, with the LD resource it
# was given, if any.
is_local_source <- function(src) inherits(src, "og_local_source")

# Whether the LD Source's LD resource is a PLINK Reference (`og_local(bfile = )`), which holds
# every pair and has no thresholds of its own.
is_plink_source <- function(ld_src) {
  inherits(opengwasr::og_ld_resource(ld_src), "og_plink_reference")
}

# A Source that names TwoSampleMR as its caller in the user agent.
twosamplemr_source <- function() {
  opengwasr::og_api(caller = paste0("TwoSampleMR/", utils::packageVersion("TwoSampleMR")))
}

alid_pattern <- "^[0-9XYMT]+:[0-9]+:[ACGT]+:[ACGT]+$"
rsid_pattern <- "^rs[0-9]+$"

# --- LD -------------------------------------------------------------------------------------------

# The Source that answers LD for a call (ADR 0006 of opengwasr: LD is chosen on the Source, and
# only there):
# - `bfile`: an LD-only Local Source on that PLINK fileset, `og_local(bfile = , plink = )`. PLINK
#   needs no API, so `src` is never read, and a call with `bfile` never makes the default Source.
# - `pop`: `src` with that LD Store, `og_api(inherit = src, ld_store = pop)`, which keeps the user's
#   URL and key. `pop` names an LD Store on the OpenGWAS API, and a Local Source answers LD only from
#   resources on this machine (its `ld_store` would be the path of a downloaded LD Store), so `pop`
#   on a Local Source is refused here (`twosamplemr_pop_local`), naming the API Source and a PLINK
#   Reference instead.
# - neither: `src` itself if it carries an LD resource (`og_api(ld_store = )`, `og_local(bfile = )`),
#   or is a Local Source, which opengwasr refuses for want of one. A plain API Source gives the
#   European LD Store for an LD matrix (`ids = NULL`), and for clumping is `src` itself, so that
#   opengwasr matches the Assigned Ancestry of the Analyses `ids` name; data whose ids name no
#   Analysis (from a file, say) have none to match, so they are refused.
ld_source <- function(src, pop, bfile, plink_bin, ids) {
  if (!is.null(bfile)) {
    return(opengwasr::og_local(bfile = bfile, plink = plink_bin))
  }
  if (identical(pop, "legacy")) {
    deprecated_argument_message(paste0(
      "`pop = \"legacy\"` is not available: the old EUR panel with a different set of markers is ",
      "gone. The European LD Store is used."
    ))
    pop <- "EUR"
  }
  if (is_local_source(src)) {
    if (!is.null(pop)) {
      pop_local_refusal(pop)
    }
    return(src)
  }
  if (is.null(pop) && !is.null(opengwasr::og_ld_resource(src))) {
    return(src)
  }
  if (!is.null(pop) || is.null(ids)) {
    return(opengwasr::og_api(inherit = src, ld_store = pop %||% "EUR"))
  }
  unknown <- tryCatch(
    {
      opengwasr::og_analyses(src, unique(as.character(ids)))
      character()
    },
    opengwasr_analysis_not_found = function(e) e$ids
  )
  if (length(unknown) == 0) {
    return(src)
  }
  stop(structure(
    class = c("twosamplemr_ld_needs_pop", "error", "condition"),
    list(
      message = paste0(
        "These ids name no OpenGWAS Analysis, so there is no Assigned Ancestry to choose an LD Store ",
        "by: ",
        paste(utils::head(unknown, 5), collapse = ", "),
        if (length(unknown) > 5) ", ...",
        ". Give `pop`: an ancestry (\"EUR\", \"AFR\", \"EAS\", \"SAS\", \"AMR\") or an LD ",
        "Store id from `opengwasr::og_ld_stores()`; or `bfile`, the path of your own PLINK reference."
      ),
      call = NULL
    )
  ))
}

# `pop` on a Local Source: it names an LD Store on the OpenGWAS API, which a Local Source never reads.
pop_local_refusal <- function(pop) {
  shown <- if (is.character(pop) && length(pop) == 1) pop else "<pop>"
  stop(structure(
    class = c("twosamplemr_pop_local", "error", "condition"),
    list(
      message = paste0(
        "`pop` names an LD Store on the OpenGWAS API, but `src` is a Local Source, which answers LD ",
        "only from resources on this machine. To clump your data with the API's LD Store, use the API ",
        "Source: `clump_data(dat, pop = \"", shown, "\", src = opengwasr::og_api())`. Or give the ",
        "Local Source a PLINK Reference, `opengwasr::og_local(gwas = , bfile = )`, and leave out `pop`."
      ),
      call = NULL,
      pop = pop,
      src_kind = "local"
    )
  ))
}

# The ALID of each row of `dat`: its `SNP` if that is one, else made from its position and alleles
# (GRCh38), else looked up from its rsID; NA where none of these names the variant.
dat_alids <- function(dat, src, type = "exposure") {
  snp <- as.character(dat$SNP)
  alid <- ifelse(grepl(alid_pattern, snp), snp, NA_character_)
  column <- function(name) dat[[paste0(name, ".", type)]]
  ea <- toupper(as.character(column("effect_allele")))
  oa <- toupper(as.character(column("other_allele")))
  chr <- column("chr")
  pos <- column("pos")
  need <- is.na(alid)
  if (any(need) && !is.null(chr) && !is.null(pos) && length(ea) > 0 && length(oa) > 0) {
    made <- data.frame(
      chromosome = as.character(chr[need]),
      position = suppressWarnings(as.numeric(pos[need])),
      effect_allele = ea[need],
      other_allele = oa[need],
      stringsAsFactors = FALSE
    )
    made <- withCallingHandlers(
      opengwasr::og_alid(made, eaf = NULL, beta = NULL, z = NULL),
      message = function(m) {
        if (
          inherits(m, "opengwasr_alid_swapped") ||
            inherits(m, "opengwasr_alid_unnamed") ||
            inherits(m, "opengwasr_alid_undeclared_reference")
        ) {
          invokeRestart("muffleMessage")
        }
      }
    )
    alid[need] <- made$alid
  }
  need <- is.na(alid) & grepl(rsid_pattern, snp, ignore.case = TRUE)
  if (any(need)) {
    found <- opengwasr::og_variants_rsid(src, unique(tolower(snp[need])))
    found <- found[!duplicated(paste(found$alid, found$rsid)), , drop = FALSE]
    for (i in which(need)) {
      candidates <- found[found$rsid == tolower(snp[i]), , drop = FALSE]
      if (!is.na(ea[i]) && !is.na(oa[i])) {
        same <- (candidates$effect_allele == ea[i] & candidates$other_allele == oa[i]) |
          (candidates$effect_allele == oa[i] & candidates$other_allele == ea[i])
        candidates <- candidates[same, , drop = FALSE]
      }
      if (nrow(candidates) == 1) alid[i] <- candidates$alid
    }
  }
  alid
}

# `dat` clumped by opengwasr with the LD Source `ld_src` (ld_source()): the rows of its index
# variants, with the clumping provenance columns.
#
# `from_bfile` says the LD Source was made from TwoSampleMR's own `bfile` argument (#28), an LD-only
# Local Source. Two things follow from that, and from nothing else:
# - variants are matched to its PLINK Reference by rsID, so no ALID is needed and nothing is asked of
#   the API: each row gets a placeholder ALID, which only keys the row, and its `SNP` is its rsID.
#   Variants the reference lacks are dropped by opengwasr, with its count message;
# - `clump_r2` and `clump_kb` default to TwoSampleMR's old 0.001 and 10000, as `bfile` always did.
# A Source's own PLINK Reference (`og_local(gwas = , bfile = )`) is matched by ALID and gets no
# defaults: opengwasr requires `r2` and `kb` for it (ADR 0006; `opengwasr_plink_settings_required`).
#
# `by_analysis` is og_clump()'s: each `id.exposure` clumped on its own (TRUE), or all the rows as one
# set (FALSE). The index rows are returned in the order of `dat`.
opengwasr_clump <- function(ld_src, dat, pval_column, p1, r2, kb, from_bfile, by_analysis) {
  if (from_bfile) {
    r2 <- r2 %||% 0.001
    kb <- kb %||% 10000
  }
  p <- as.numeric(dat[[pval_column]])
  p[is.na(p)] <- 1
  if (from_bfile) {
    alid <- paste0("1:", seq_len(nrow(dat)), ":A:C")
  } else {
    alid <- dat_alids(dat, ld_src)
  }
  usable <- !is.na(alid) & p <= p1
  if (!from_bfile && any(is.na(alid))) {
    message(
      sum(is.na(alid)),
      " variant(s) could not be given an ALID (they need a `SNP` that is an rsID or an ALID, ",
      "and GRCh38 `chr` and `pos`), so they are left out of the clumping: ",
      paste(utils::head(dat$SNP[is.na(alid)], 5), collapse = ", "),
      if (sum(is.na(alid)) > 5) ", ..."
    )
  }
  d <- data.frame(
    alid = alid[usable],
    z = abs(stats::qnorm(pmax(p[usable], 1e-300) / 2)),
    id = as.character(dat$id.exposure[usable]),
    row = which(usable),
    stringsAsFactors = FALSE
  )
  if (from_bfile) {
    rsid <- as.character(dat$SNP[usable])
    d$rsid <- ifelse(grepl(rsid_pattern, rsid, ignore.case = TRUE), rsid, NA_character_)
  }
  out <- opengwasr::og_clump(ld_src, d, r2 = r2, kb = kb, by_analysis = by_analysis)
  index <- sort(unique(out$row))
  clumped <- dat[index, , drop = FALSE]
  rownames(clumped) <- NULL
  for (column in c("clump_ld", "clump_ld_choice", "clump_r2", "clump_kb")) {
    clumped[[paste0(column, ".exposure")]] <- out[[column]][match(index, out$row)]
  }
  clumped
}

# The LD matrix of `snps` from the LD Source `ld_src` (ld_source()), as a plain matrix.
opengwasr_ld_matrix <- function(ld_src, snps, fill_below_t_far) {
  m <- opengwasr::og_ld_matrix(ld_src, snps, fill_below_t_far = fill_below_t_far)
  # A plain matrix, as before: opengwasr's provenance and class are left behind.
  m <- unclass(m)
  attr(m, "og_provenance") <- NULL
  m
}

# --- Helpers --------------------------------------------------------------------------------------

# The columns of an opengwasr result that TwoSampleMR's formats have no slot for, copied across by
# Analysis and variant after format_data() has dropped them. `from` is the opengwasr data frame,
# `to` the formatted one, and `index` says which row of `from` each row of `to` came from.
carry_columns <- function(to, from, index, columns, type) {
  for (column in columns) {
    if (column %in% names(from)) {
      to[[paste0(column, ".", type)]] <- from[[column]][index]
    }
  }
  to
}

# opengwasr names variants `rsid`, which is NA for a variant without an rsID; `alid` always exists.
variant_name <- function(d, requested = NULL) {
  name <- ifelse(is.na(d$rsid), d$alid, d$rsid)
  if (!is.null(requested)) {
    # Someone who asked for an ALID gets that ALID back, even where the Analysis also knows an rsID.
    name[d$alid %in% requested] <- d$alid[d$alid %in% requested]
  }
  as.character(name)
}
