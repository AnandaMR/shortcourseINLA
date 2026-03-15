## remotes::install_github("rfsaldanha/microdatasus")

library(sf)
library(microdatasus)
library(INLA)

if(file.exists("data/sinasc.rds")) {
    sinasc <- readRDS("data/sinasc.rds")
} else {
    sinasc <- fetch_datasus(
        year_start = 2024, 
        year_end = 2024, 
        uf = "MG", 
        information_system = "SINASC"
    )
    saveRDS(sinasc, "data/sinasc.rds")
}

dim(sinasc)

names(sinasc)

## MAPA
##
if(file.exists("maps/mg_mun.shp")) {
    mg_mun <- st_read("maps/mg_mun.shp")
} else {
    library(geobr)
    mg_mun <- read_municipality(
        code_muni = "MG", 
        year = 2020, 
        showProgress = TRUE
    )
    st_write(mg_mun, "maps/mg_mun.shp")
}

if(!file.exists("grafo")) {
    library(spdep)
    viz <- poly2nb(mg_mun)
    nb2INLA(file = "grafo", nb = viz) 
}

table(sinasc$APGAR1)

dat0 <- data.frame(
    ap1 = (as.integer(sinasc$APGAR1)<=7)+0L,
    ap5 = (as.integer(sinasc$APGAR5)<=7)+0L,
    peso= as.numeric(sinasc$PESO)/1000
)

dat0$area <- pmatch(
    x = sinasc$CODMUNRES,
    table = substr(mg_mun$code_mn, 1,6),
    duplicates.ok = TRUE)

summary(dat0)

hist(dat0$peso)

dat0$pesoCat = findInterval(dat0$peso, c(-Inf, 2:9/2, Inf))
table(dat0$pesoCat)

## eta_i = f_k(peso_i) + u_i

fs1 <- ap1 ~ 0 + 
    f(pesoCat, model = 'rw2', constr = FALSE,
      hyper = list(theta = list('loggamma', param = c(2, 2)))) +
    f(area, model = 'iid', constr = TRUE,
      hyper = list(theta = list('loggamma', param = c(2, 2))))

fit1 <- inla(
    formula = fs1, 
    family  = 'binomial',
    data    = dat0, 
    control.predictor = list(link = 1)
)

fit1$summary.fixed

inla.priors.used(fit1)

head(INLA:::plot.inla,3)

plot(fit1, T, F, T, T, F, F, F, plot.prior = TRUE)

## PC-prec prior: P(sigma > sigma_0) = p
pcprec <- list(prior = 'pc.prec', param = c(0.5, 0.01))

fs_pc <- ap1 ~ 0 +
    f(pesoCat, model = 'rw2', constr = FALSE,
      scale.model = TRUE,
      hyper = list(prec = pcprec)) +
    f(area, model = 'iid', hyper = list(prec = pcprec))

fit_pc <- inla(
    formula = fs_pc, 
    family  = 'binomial',
    data    = dat0, 
    control.predictor = list(link = 1)
)

plot(fit_pc, T, F, T, T, F, F, F, plot.prior = TRUE)

n <- nrow(dat0)
str(dat0)

stk1 <- inla.stack(
    tag = 'ap1',
    data = list(ap1 = dat0$ap1, link_id = rep(1, n)),
    effect = list(
        data.frame(pesoCat1 = dat0$pesoCat,
                   area1 = dat0$area)),
    A = list(1)
)

stk5 <- inla.stack(
    tag = 'ap5',
    data = list(ap5 = dat0$ap5, link_id = rep(2, n)),
    effect = list(
        data.frame(pesoCat5 = dat0$pesoCat,
                   area5 = dat0$area)),
    A = list(1)
)

## joint data stack
stkj <- inla.stack(stk1, stk5)

ffj <- list(ap1, ap5) ~ 0 +
    f(pesoCat1,  model = 'rw2', constr = FALSE,
      scale.model = TRUE, hyper = list(prec = pcprec)) +
    f(area1, model = 'besagproper', graph = "grafo",
      hyper = list(prec = pcprec)) +
    f(pesoCat5,  model = 'rw2', constr = FALSE,
      scale.model = TRUE, hyper = list(prec = pcprec)) +
    f(area5, copy = 'area1', fixed = FALSE)


fitj <- inla(
    formula = ffj, family = rep("binomial", 2),
    data = inla.stack.data(stkj),
    control.predictor = list(
        A = inla.stack.A(stkj),
        link = link_id)
)

fitj$cpu.used

fitj$summary.hy

stk1r <- inla.stack(
    tag = 'ap1',
    data = list(ap1 = dat0$ap1, link_id = rep(1, n)),
    effect = list(
        data.frame(pesoCat = dat0$pesoCat, r = rep(1, n),
                   area1 = dat0$area)),
    A = list(1)
)

stk5r <- inla.stack(
    tag = 'ap5',
    data = list(ap5 = dat0$ap5, link_id = rep(2, n)),
    effect = list(
        data.frame(pesoCat = dat0$pesoCat, r = rep(2, n),
                   area5 = dat0$area)),
    A = list(1)
)

stkjr <- inla.stack(stk1r, stk5r)

ffjr <- list(ap1, ap5) ~ 0 +
    f(pesoCat,  model = 'rw2', replicate = r, constr = FALSE,
      scale.model = TRUE, hyper = list(prec = pcprec)) +
    f(area1, model = 'besagproper', graph = "grafo",
      hyper = list(prec = pcprec)) +
    f(area5, copy = 'area1', fixed = FALSE)


fitjr <- inla(
    formula = ffjr, family = rep("binomial", 2),
    data = inla.stack.data(stkjr),
    control.predictor = list(
        A = inla.stack.A(stkjr),
        link = link_id)
)

fitjr$cpu.used

fitjr$summary.hy

plot(fitjr, T, F, T, T, F, F, F, plot.prior = TRUE)

inla.priors.used(fitjr)
