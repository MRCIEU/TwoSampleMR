docs:
    Rscript -e "devtools::document()"
check: docs
    Rscript -e "devtools::check()"
install: docs
    Rscript -e "devtools::install(build_vignettes = TRUE)"
dev:
    Rscript -e "pak::local_install_dev_deps()"
