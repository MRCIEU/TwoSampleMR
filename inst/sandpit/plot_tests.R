library(devtools)
load_all()

ed <- extract_instruments(2)
ed1 <- ed[1,]
ed2 <- ed[1:2,]
ed3 <- ed[1:3,]

od <- extract_outcome_data(ed$SNP, 7)
od1 <- extract_outcome_data(ed1$SNP, 7)
od2 <- extract_outcome_data(ed2$SNP, 7)
od3 <- extract_outcome_data(ed3$SNP, 7)

dat <- harmonise_data(ed, od)
dat1 <- harmonise_data(ed1, od1)
dat2 <- harmonise_data(ed2, od2)
dat3 <- harmonise_data(ed3, od3)

mr_heterogeneity(dat)
mr_heterogeneity(dat1)
mr_heterogeneity(dat2)
mr_heterogeneity(dat3)

mrs <- mr_singlesnp(dat)
mrs1 <- mr_singlesnp(dat1)
mrs2 <- mr_singlesnp(dat2)
mrs3 <- mr_singlesnp(dat3)

mr_forest_plot(mrs)
mr_forest_plot(mrs1)
mr_forest_plot(mrs2)
mr_forest_plot(mrs3)

mr_funnel_plot(mrs)
mr_funnel_plot(mrs1)
mr_funnel_plot(mrs2)
mr_funnel_plot(mrs3)

mr_scatter_plot(mr(dat), dat)
mr_scatter_plot(mr(dat1), dat1)
mr_scatter_plot(mr(dat2), dat2)
mr_scatter_plot(mr(dat3), dat3)


mrl <- mr_leaveoneout(dat)
mrl1 <- mr_leaveoneout(dat1)
mrl2 <- mr_leaveoneout(dat2)
mrl3 <- mr_leaveoneout(dat3)

mr_leaveoneout_plot(mrl)
mr_leaveoneout_plot(mrl1)
mr_leaveoneout_plot(mrl2)
mr_leaveoneout_plot(mrl3)

