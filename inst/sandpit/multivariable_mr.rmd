
```{r}
library(devtools)
library(knitr)
load_all()


lipids <- mv_extract_exposures(c(299,300,302))
chd <- extract_outcome_data(lipids$SNP, 7)
control <- mv_harmonise_data(lipids, chd)
```


```{r}
kable(mv_residual(control, intercept=TRUE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_residual(control, intercept=FALSE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_residual(control, intercept=TRUE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_residual(control, intercept=FALSE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_multiple(control, intercept=TRUE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_multiple(control, intercept=FALSE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_multiple(control, intercept=TRUE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_multiple(control, intercept=FALSE, instrument_specific=FALSE)$result)
```


```{r}
a <- mv_extract_exposures(c("UKB-a:196", 1001))
b <- extract_outcome_data(a$SNP, 297)
dat <- mv_harmonise_data(a, b)
```

```{r}
kable(mv_residual(dat, intercept=TRUE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_residual(dat, intercept=FALSE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_residual(dat, intercept=TRUE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_residual(dat, intercept=FALSE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_multiple(dat, intercept=TRUE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_multiple(dat, intercept=FALSE, instrument_specific=TRUE)$result)
```

```{r}
kable(mv_multiple(dat, intercept=TRUE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_multiple(dat, intercept=FALSE, instrument_specific=FALSE)$result)
```

```{r}
kable(mv_ivw(dat)$result)
```

```{r}
kable(mv_basic(dat)$result)
```


