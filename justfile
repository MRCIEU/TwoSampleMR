docs:
    Rscript -e "devtools::document()"
check: docs
    Rscript -e "devtools::check()"
test:
    Rscript -e "devtools::test()"
test-server:
    TWOSAMPLEMR_ENABLE_OPENGWAS_TESTS=TRUE Rscript -e "devtools::test()"
install: docs
    Rscript -e "devtools::install(build_vignettes = TRUE)"
dev:
    Rscript -e "pak::local_install_dev_deps()"
readme:
    Rscript -e "rmarkdown::render('README.Rmd', output_options = list(html_preview = FALSE))"
