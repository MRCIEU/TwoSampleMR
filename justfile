docs:
    R -e "devtools::document()"
check: docs
    R -e "devtools::check()"
install: docs
    R -e "devtools::install(build_vignettes = TRUE)"
dev:
    R -e "pak::local_install_dev_deps()"
