skip_if_opengwas_tests_disabled <- function() {
  if (!identical(Sys.getenv("TWOSAMPLEMR_ENABLE_OPENGWAS_TESTS"), "TRUE")) {
    skip("Tests requiring OpenGWAS.")
  } else {
    invisible()
  }
}
