Ghuber <- function (u, k = 30, deriv = 0)
{
    if (!deriv)
    {
        return(pmin(1, k/abs(u)))
    } else {
        return(abs(u) <= k)
    }
}

#' Univariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param Z summary Z-statistics for M variants
#' @param r2 average reference LD scores for M variants
#' @param N GWAS sample size for each variant (could be different across variants)
#' @param W variant weight
#'
#' @keywords internal
#' @return model fit
ldsc_h2_internal <- function(Z, r2, N, W=NULL)
{
    if (is.null(W))
    {
        W <- rep(1, length(Z))
    }
    tau <- (mean(Z^2) - 1) / mean(N * r2)
    Wv <- 1 / (1 + tau * N * r2)^2
    id <- which(Z^2 > 30)
    if (length(id) > 0)
    {
        Wv[id] <- sqrt(Wv[id])
    }
    mod <- MASS::rlm(I(Z^2) ~ I(N * r2), weight = W * Wv,
        psi = Ghuber, k = 30)
    return(summary(mod))
}


#' Bivariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param Zs Mx2 matrix of summary Z-statistics for M variants from two GWAS
#' @param r2 average reference LD scores for M variants
#' @param N1 sample size for the 1st GWAS
#' @param N2 sample size for the 2nd GWAS
#' @param Nc overlapped sample size between the two GWAS
#' @param W variant weight
#' @param h1 hsq for trait 1
#' @param h2 hsq for trait 2
#'
#' @return List of models
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations across human diseases and traits. Nat. Genet. 47, 1236–1241.
#'
#' Guo,B. and Wu,B. (2018) Principal component based adaptive association test of multiple traits using GWAS summary statistics. bioRxiv 269597; doi: 10.1101/269597
#'
#' Gua,B. and Wu,B. (2019) Integrate multiple traits to detect novel trait-gene association using GWAS summary data with an adaptive test approach. Bioinformatics. 2019 Jul 1;35(13):2251-2257. doi: 10.1093/bioinformatics/bty961. 
#'
#' https://github.com/baolinwu/MTAR 
#' @keywords internal
ldsc_rg_internal <- function(Zs, r2, h1, h2, N1, N2, Nc=0, W=NULL)
{
    if(is.null(W))
    {
        W = rep(1,length(r2))
    }

    Y <- Zs[,1] * Zs[,2]

    X <- (sqrt(N1) * sqrt(N2) + sqrt(Nc/N1 * N2)) * r2
    N1r2 <- N1 * r2
    N2r2 <- N2 * r2
    r0 <- 0

    ## 1st round
    if(any(Nc > 0))
    {
        rcf <- as.vector(MASS::rlm(Y ~ X, psi = Ghuber)$coef)
        r0 <- rcf[1]
        gv <- rcf[-1]
    } else {
        gv <- as.vector(MASS::rlm(Y ~ X-1, psi = Ghuber)$coef)
    }

    ## 2nd round
    Wv <- 1 / ((h1 * N1r2 + 1) * (h2 * N2r2 + 1) + (X * gv + r0)^2)
    id <- which(abs(Zs[,1] * Zs[,2]) > 30)
    if(length(id) > 0)
    {
        Wv[id] <- sqrt(Wv[id])
    }

    if(any(Nc > 0))
    {
        rcf <- MASS::rlm(Y ~ X, weight = W * Wv, psi = Ghuber, k = 30)
    } else {
        rcf <- MASS::rlm(Y ~ X - 1, weight = W * Wv, psi = Ghuber, k = 30)
    }
    return(summary(rcf))
}


#' Univariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param id ID to analyse
#' @param ancestry ancestry of traits 1 and 2 (AFR, AMR, EAS, EUR, SAS) or 'infer' (default) in which case it will try to guess based on allele frequencies
#' @param snpinfo Output from ieugwasr::afl2_list("hapmap3"), or NULL for it to be done automatically
#' @param splitsize How many SNPs to extract at one time. Default=20000
#'
#' @export
#' @return model fit
#' @references
#' Bulik-Sullivan,B.K. et al. (2015) An atlas of genetic correlations across human diseases and traits. Nat. Genet. 47, 1236–1241.
#'
#' Guo,B. and Wu,B. (2018) Principal component based adaptive association test of multiple traits using GWAS summary statistics. bioRxiv 269597; doi: 10.1101/269597
#'
#' Gua,B. and Wu,B. (2019) Integrate multiple traits to detect novel trait-gene association using GWAS summary data with an adaptive test approach. Bioinformatics. 2019 Jul 1;35(13):2251-2257. doi: 10.1093/bioinformatics/bty961. 
#'
#' https://github.com/baolinwu/MTAR 
ldsc_h2 <- function(id, ancestry="infer", snpinfo = NULL, splitsize=20000)
{
    if(is.null(snpinfo))
    {
        snpinfo <- ieugwasr::afl2_list("hapmap3")
    }

    snpinfo <- snpinfo %>%
        dplyr::filter(complete.cases(.))

    d <- extract_split(snpinfo$rsid, id, splitsize) %>%
        ieugwasr::fill_n() %>%
        dplyr::mutate(z = beta / se) %>%
        dplyr::select(rsid, z = z, n = n, eaf) %>%
        dplyr::filter(complete.cases(.))

    stopifnot(nrow(d) > 0)

    if(ancestry == "infer")
    {
        ancestry <- ieugwasr::infer_ancestry(d, snpinfo)$pop[1]
    }

    d <- snpinfo %>% 
        dplyr::select(rsid, l2=paste0("L2.", ancestry)) %>%
        dplyr::inner_join(., d, by="rsid") %>%
        dplyr::filter(complete.cases(.))

    return(ldsc_h2_internal(d$z, d$l2, d$n))
}

#' Bivariate LDSC
#'
#' Imported here to help estimate sample overlap between studies
#'
#' @param id1 ID 1 to analyse
#' @param id2 ID 2 to analyse
#' @param ancestry ancestry of traits 1 and 2 (AFR, AMR, EAS, EUR, SAS) or 'infer' (default) in which case it will try to guess based on allele frequencies
#' @param snpinfo Output from ieugwasr::afl2_list("hapmap3"), or NULL for it to be done automatically
#' @param splitsize How many SNPs to extract at one time. Default=20000
#'
#' @export
#' @return model fit
ldsc_rg <- function(id1, id2, ancestry="infer", snpinfo = NULL, splitsize=20000)
{
    if(is.null(snpinfo))
    {
        snpinfo <- ieugwasr::afl2_list("hapmap3")
    }

    x <- extract_split(snpinfo$rsid, c(id1, id2), splitsize)
    d1 <- subset(x, id == id1) %>%
        ieugwasr::fill_n() %>%
        dplyr::mutate(z = beta / se) %>%
        dplyr::select(rsid, z1 = z, n1 = n, eaf) %>%
        dplyr::filter(complete.cases(.))

    stopifnot(nrow(d1) > 0)

    d2 <- subset(x, id == id2) %>%
        ieugwasr::fill_n() %>%
        dplyr::mutate(z = beta / se) %>%
        dplyr::select(rsid, z2 = z, n2 = n, eaf) %>%
        dplyr::filter(complete.cases(.))

    stopifnot(nrow(d2) > 0)

    if(ancestry == "infer")
    {
        ancestry1 <- ieugwasr::infer_ancestry(d1, snpinfo)
        ancestry2 <- ieugwasr::infer_ancestry(d2, snpinfo)
        if(ancestry1$pop[1] != ancestry2$pop[1])
        {
            stop("d1 ancestry is ", ancestry1$pop[1], " and d2 ancestry is ", ancestry2$pop[1])
        }
        ancestry <- ancestry1$pop[1]
    }

    d1 <- snpinfo %>% 
        dplyr::select(rsid, l2=paste0("L2.", ancestry)) %>%
        dplyr::inner_join(., d1, by="rsid")

    d2 <- snpinfo %>% 
        dplyr::select(rsid, l2=paste0("L2.", ancestry)) %>%
        dplyr::inner_join(., d2, by="rsid")

    h1 <- ldsc_h2_internal(d1$z1, d1$l2, d1$n1, d1$l2)
    h2 <- ldsc_h2_internal(d2$z2, d2$l2, d2$n2, d1$l2)

    dat <- dplyr::inner_join(d1, d2, by="rsid") %>%
        dplyr::mutate(
            l2 = l2.x,
            n1 = as.numeric(n1),
            n2 = as.numeric(n2),
            rhs = l2 * sqrt(n1 * n2)
        )

    gcov <- dat %>%
        {
            ldsc_rg_internal(
                Zs = cbind(.$z1, .$z2),
                r2 = .$l2,
                h1 = h1$coefficients[2,1] * nrow(d1),
                h2 = h2$coefficients[2,1] * nrow(d2),
                N1 = .$n1,
                N2 = .$n2,
                W = .$l2
            )
        }
    return(list(
        gcov = gcov,
        h1=h1,
        h2=h2,
        rg = (gcov$coefficients[1,1] * nrow(dat)) / sqrt(h1$coefficients[2,1] * nrow(d1) * h2$coefficients[2,1] * nrow(d2))
    ))
}


extract_split <- function(snplist, id, splitsize=20000)
{
    nsplit <- round(length(snplist)/splitsize)
    split(snplist, 1:nsplit) %>%
        pbapply::pblapply(., function(x)
        {
            ieugwasr::associations(x, id, proxies=FALSE)
        }) %>% dplyr::bind_rows()
}
