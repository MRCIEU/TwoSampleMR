skip_if_opengwas_tests_disabled <- function() {
  if (!identical(Sys.getenv("TWOSAMPLEMR_ENABLE_OPENGWAS_TESTS"), "TRUE")) {
    skip("Tests requiring OpenGWAS.")
  } else {
    invisible()
  }
}

# opengwasr answers TwoSampleMR's data calls. Its fake Data API is internal to opengwasr's own tests,
# so these tests replace opengwasr's functions at the boundary TwoSampleMR calls, with small data
# frames recorded from real opengwasr output (see fixtures/record-fixtures.R).
fixture <- function(name) readRDS(test_path("fixtures", paste0(name, ".rds")))

bmi_id <- "OGS-00009:ukb-b-19953"
cad_id <- "OGS-00011:GCST003116"

# The calls opengwasr received, most recent last, are in `$calls`.
local_fake_opengwasr <- function(empty = character(), env = parent.frame()) {
  analyses <- fixture("analyses")
  tophits <- list(
    "OGS-00009:ukb-b-19953" = fixture("tophits-bmi"),
    "OGS-00011:GCST003116" = fixture("tophits-cad")
  )
  canonical <- c(
    "ukb-b-19953" = bmi_id,
    "ebi-a-GCST003116" = cad_id,
    stats::setNames(c(bmi_id, cad_id), c(bmi_id, cad_id))
  )
  associations <- fixture("associations-cad")
  state <- new.env()
  state$calls <- list()
  note <- function(fn, args) state$calls[[length(state$calls) + 1]] <- list(fn = fn, args = args)

  testthat::local_mocked_bindings(
    og_api = function(...) structure(list(), class = "fake_source"),
    og_tophits = function(
      src,
      analyses,
      threshold = 5e-8,
      limit = 1000,
      observed_only = FALSE,
      clump = TRUE,
      ld = NULL,
      r2 = NULL,
      kb = NULL,
      original_units = FALSE,
      releases = NULL
    ) {
      note(
        "og_tophits",
        list(
          src = src,
          analyses = analyses,
          threshold = threshold,
          limit = limit,
          clump = clump,
          ld = ld,
          r2 = r2,
          kb = kb
        )
      )
      rows <- lapply(analyses, function(a) {
        id <- canonical[[a]]
        if (id %in% empty) {
          return(NULL)
        }
        out <- tophits[[id]]
        if (!clump) {
          out <- out[, !grepl("^clump_", names(out))]
        }
        out
      })
      out <- do.call(rbind, rows)
      if (is.null(out)) tophits[[bmi_id]][0, ] else out
    },
    og_associations = function(
      src,
      variants,
      analyses,
      observed_only = FALSE,
      original_units = FALSE,
      releases = NULL
    ) {
      note("og_associations", list(src = src, variants = variants, analyses = analyses))
      # BMI's Associations of its own instruments are its recorded top hits.
      clump_cols <- c("clump_ld", "clump_ld_choice", "clump_r2", "clump_kb", "clump_size")
      bmi <- tophits[[bmi_id]][, setdiff(names(tophits[[bmi_id]]), clump_cols)]
      served <- list(bmi, associations)
      names(served) <- c(bmi_id, cad_id)
      rows <- lapply(analyses, function(a) {
        x <- served[[canonical[[a]]]]
        x[x$rsid %in% variants | x$alid %in% variants, ]
      })
      do.call(rbind, rows)
    },
    og_analyses = function(src, ids = NULL, releases = NULL) {
      note("og_analyses", list(src = src, ids = ids))
      if (is.null(ids)) {
        return(analyses)
      }
      row <- match(ids, analyses$requested_id)
      row[is.na(row)] <- match(ids, analyses$id)[is.na(row)]
      unknown <- ids[is.na(row)]
      found <- analyses[row[!is.na(row)], ]
      found$requested_id <- ids[!is.na(row)]
      if (length(unknown) > 0) {
        stop(structure(
          class = c("opengwasr_analysis_not_found", "error", "condition"),
          list(message = "Identifiers that name no Analysis", call = NULL, ids = unknown)
        ))
      }
      found
    },
    .package = "opengwasr",
    .env = env
  )
  state
}

calls_to <- function(state, fn) Filter(function(x) x$fn == fn, state$calls)

# What og_clump() answers: the rows of `data` named in `keep` (by row of the data it is given, so
# the default keeps the first row of every window of 1 Mb), with its clumping columns; and what
# og_ld_matrix() answers: an og_ld_matrix named by ALID.
local_fake_ld <- function(keep = NULL, env = parent.frame()) {
  state <- new.env()
  state$calls <- list()
  testthat::local_mocked_bindings(
    og_clump = function(src, data, ld = NULL, r2 = NULL, kb = NULL) {
      state$calls[[length(state$calls) + 1]] <- list(
        fn = "og_clump",
        src = src,
        data = data,
        ld = ld,
        r2 = r2,
        kb = kb
      )
      rows <- if (identical(keep, "all")) {
        seq_len(nrow(data))
      } else if (is.null(keep)) {
        seq(1, nrow(data), by = 2)
      } else {
        keep
      }
      out <- data[rows, , drop = FALSE]
      out$clump_ld <- "hgdp1kgp-v2-EUR-k10"
      out$clump_ld_choice <- "matched"
      out$clump_r2 <- r2 %||% 0.0319
      out$clump_kb <- kb %||% 10000
      out$clump_size <- 2L
      out
    },
    og_ld_matrix = function(src, variants, ld = NULL, fill_below_t_far = FALSE) {
      state$calls[[length(state$calls) + 1]] <- list(
        fn = "og_ld_matrix",
        # TwoSampleMR's default Source needs credentials, which these tests don't have.
        src = tryCatch(src, opengwasr_auth = function(e) NULL),
        variants = variants,
        ld = ld,
        fill = fill_below_t_far
      )
      m <- diag(length(variants))
      dimnames(m) <- list(variants, variants)
      structure(
        m,
        class = c("og_ld_matrix", "matrix", "array"),
        og_provenance = list(ld = "hgdp1kgp-v2-EUR-k10")
      )
    },
    .package = "opengwasr",
    .env = env
  )
  state
}

# The ALID of each instrument, from the recorded top hits (an ALID is its own).
alid_of <- function(snp) {
  hits <- fixture("tophits-bmi")
  ifelse(snp %in% hits$rsid, hits$alid[match(snp, hits$rsid)], snp)
}

# #28: nothing in the PLINK path calls the API, so it works with no credentials and no Source.
local_no_credentials <- function(env = parent.frame()) {
  withr::local_envvar(
    OPENGWAS_URL = NA,
    OPENGWAS_DEV_KEY = NA,
    OPENGWAS_DEV_USER = NA,
    OPENGWAS_JWT = NA,
    .local_envir = env
  )
  withr::local_options(opengwasr.url = NULL, .local_envir = env)
}

# --- Local Source ---------------------------------------------------------------------------------

# opengwasr's example Store Releases (one per GWAS Store: `hyb` holds two, so its `v1` is named) and
# its synthetic PLINK Reference, whose rsIDs include some of their top hits.
example_releases <- function() {
  stores <- system.file("stores", package = "opengwasr")
  c(file.path(stores, c("dense", "ragged")), file.path(stores, "hyb", "v1"))
}
example_bfile <- function() file.path(system.file("plink", package = "opengwasr"), "ref")

# As opengwasr's own: the local tests need a Python with ogquery named by RETICULATE_PYTHON, and
# never set Python up themselves.
skip_if_no_local <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  has_local <- "og_local" %in% getNamespaceExports("opengwasr")
  if (!has_local || !nzchar(system.file("stores", package = "opengwasr"))) {
    skip("This opengwasr has no Local Source, or no example Store Releases.")
  }
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    skip("RETICULATE_PYTHON isn't set; local tests never set Python up themselves.")
  }
  if (is.null(example_local())) skip("The Python in RETICULATE_PYTHON has no ogquery.")
  invisible()
}

skip_if_no_local_plink <- function() {
  skip_if_no_local()
  if (is.null(example_local(plink = TRUE))) {
    skip("No PLINK binary found (genetics.binaRies, or plink on the PATH).")
  }
  invisible()
}

# A Local Source over the example Store Releases, with the example PLINK Reference as its LD resource
# (`plink = TRUE`) or with none, opened once per test process; NULL where it can't be opened (no
# ogquery, or no PLINK binary).
example_local <- local({
  opened <- list()
  function(plink = FALSE) {
    key <- if (plink) "plink" else "none"
    if (!key %in% names(opened)) {
      opened[key] <<- list(tryCatch(
        if (plink) {
          opengwasr::og_local(example_releases(), bfile = example_bfile())
        } else {
          opengwasr::og_local(example_releases())
        },
        opengwasr_local_unavailable = function(e) NULL,
        opengwasr_plink_not_found = function(e) NULL
      ))
    }
    opened[[key]]
  }
})

# --- The Source's LD resource (ADR 0006) ------------------------------------------------------------

# A real API Source, made with placeholder credentials: making one calls nothing. Make it before
# local_fake_opengwasr(), which replaces og_api().
api_source <- function(ld_store = NULL) {
  withr::with_envvar(
    c(OPENGWAS_URL = "https://api.opengwas.invalid", OPENGWAS_DEV_KEY = "test-key"),
    opengwasr::og_api(ld_store = ld_store)
  )
}

# The LD that opengwasr uses for a call given `src` and `ld`: the call's `ld`, else the Source's LD
# resource, else (NULL) each Analysis's Assigned Ancestry. Asked of opengwasr's own internal
# ld_resource(), so the tests check what opengwasr would do with what TwoSampleMR sends.
ld_used <- function(src, ld) {
  utils::getFromNamespace("ld_resource", "opengwasr")(src, ld)$ld
}
