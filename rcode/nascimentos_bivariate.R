library(INLA)
library(sf)

if(file.exists("data/sinasc.rds")) {
    sinasc <- readRDS("data/sinasc.rds")
} else {
## remotes::install_github("rfsaldanha/microdatasus")
    library(microdatasus)
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

dim(sinasc)

names(sinasc)

table(sinasc$CONSULTAS)

table(sinasc$APGAR1)
table(as.numeric(sinasc$APGAR1)<=6)

summary(as.numeric(sinasc$PESO)/1000)

table(sinasc$APGAR1,
      sinasc$APGAR5)

dat0 <- data.frame(
    ap1 = (as.integer(sinasc$APGAR1)<=7)+0L,
    ap5 = (as.integer(sinasc$APGAR5)<=7)+0L,
    peso= as.numeric(sinasc$PESO)/1000,
    area= pmatch(
        x = sinasc$CODMUNRES,
        table = substr(mg_mun$code_mn, 1,6),
        duplicates.ok = TRUE)
)

summary(dat0)

fit0a1 <- inla(
    formula = ap1 ~ peso,
    family = "binomial",
    data = dat0,
    control.predictor = list(link = 1))

fit0a5 <- inla(
    formula = ap5 ~ peso,
    family = "binomial",
    data = dat0,
    control.predictor = list(link = 1))

fit0a1$summary.fixed
fit0a5$summary.fixed

## m1: eta_i = beta_0 + beta_1 * peso_i + u_mun(i)
##     P(y_i=1) = p_i = 1/(1+exp(-eta_i))
fit1a1 <- inla(
    formula = ap1 ~ peso + f(area, model = 'iid'),
    family = "binomial",
    data = dat0,
    control.predictor = list(link = 1))
fit1a5 <- inla(
    formula = ap5 ~ peso + f(area, model = 'iid'),
    family = "binomial",
    data = dat0,
    control.predictor = list(link = 1))

fit1a1$summary.fixed
fit1a5$summary.fixed

fit1a1$summary.hyperpar
fit1a5$summary.hyperpar

n <- nrow(dat0)
dat2 <- data.frame(
    ap1 = c(dat0$ap1, rep(NA, n)),
    ap5 = c(rep(NA, n), dat0$ap5),
    b0_a1 = c(rep(1., n), rep(NA, n)), ##rep(1:0, c(n,n)),
    b0_a5 = c(rep(NA, n), rep(1., n)), ##rep(0:1, c(n,n)),
    peso_a1 = c(dat0$peso, rep(NA, n)),
    peso_a5 = c(rep(NA, n), dat0$peso)
)

fit0_b <- inla(
    formula = list(ap1, ap5) ~ -1 +
        b0_a1 + b0_a5 + peso_a1 + peso_a5,
    family = rep("binomial", 2),
    data = dat2,
    control.predictor = list(link = 1))

fit0a1$summary.fixed
fit0a5$summary.fixed

fit0_b$summary.fixed

dat2$area1 <- c(dat0$area, rep(NA,n))
dat2$area2 <- c(rep(NA,n), dat0$area)

fit1_b <- inla(
    formula = list(ap1, ap5) ~ -1 +
        b0_a1 + b0_a5 + peso_a1 + peso_a5 +
        f(area1, model = "iid") +
        f(area2, model = "iid"),
    family = rep("binomial", 2),
    data = dat2,
    control.predictor = list(link = 1))

rbind(fit1a1$cpu.used, fit1a5$cpu.used, fit1_b$cpu.used)

fit1a1$summary.fixed
fit1a5$summary.fixed

fit1_b$summary.fixed

## eta1_i: a_1 + beta_1 * peso_i + u_mun(i)
##   u_mun(i) ~ iid
## eta2_i: a_2 + beta_2 * peso_i + \beta_u * u_mun(i)
##   NOTA: \beta_u coef de associacao

ff2b <- list(ap1, ap5) ~ -1 +
        b0_a1 + b0_a5 + peso_a1 + peso_a5 +
        f(area1, model = "iid") +
        f(area2, copy = 'area1', fixed = FALSE)

fit2_b <- inla(
    formula = ff2b,
    family = rep("binomial", 2),
    data = dat2,
    control.predictor = list(link = 1)
)

fit2_b$summary.hyperpar

fit3_b <- inla(
    formula = list(ap1, ap5) ~ -1 +
        b0_a1 + b0_a5 + peso_a1 + peso_a5 +
        f(area1, model = "besagproper", graph = "grafo") +
        f(area2, copy = 'area1', fixed = FALSE),
    family = rep("binomial", 2),
    data = dat2,
    control.predictor = list(link = 1)
)

fit3_b$summary.fixed

fit3_b$summary.hyperpar

## copiar todo o preditor linear eta1 em eta2
## eta1_i: a_1 + beta_1 * peso_i + u_mun(i)
## defina
##   v_i = beta_1 * peso_i + u_mun(i)
##   0.0 = beta_1 * peso_i + u_mun(i) - v_i
## eta2_i: a_2 + beta_2 * peso_i + beta_v * v_i
##   NOTA: \beta_v coef de associacao

data3 <- data.frame(
    ap1 = c(dat0$ap1, rep(NA, n),  rep(NA, n)),
    ap5 = c(rep(NA, n), dat0$ap5,  rep(NA, n)),
    zero= c(rep(NA, n), rep(NA,n),  rep(0, n)),
    lnk = rep(1:3, each = n),
    b0_a1 = rep(c(1,0,0), c(n,n,n)),
    b0_a5 = rep(c(0,1,0), c(n,n,n)),
    peso1 = c(dat0$peso, rep(NA, n), dat0$peso),
    peso5 = c(rep(NA, n), dat0$peso, rep(NA, n)), ## eff. differencial de peso em ap5
    area  = c(dat0$area,  rep(NA, n), dat0$area),
    vi    = c(rep(NA, n), rep(NA, n), 1:n),
    vi_w  = rep(c(0,-1), c(2*n, n)),
    vi_cp = c(rep(NA, n), 1:n, rep(NA, n))
)

fff <- list(ap1, ap5, zero) ~ -1 +
    b0_a1 + b0_a5 + peso1 + peso5 +
    f(area, model = "besagproper", graph = "grafo") +
    f(vi, vi_w, model = "iid",
      hyper = list(theta = list(initial = -20, fixed = TRUE))) +
    f(vi_cp, copy = 'vi', fixed = FALSE)


fit3z <- inla(
    formula = fff, 
    family = c(rep("binomial", 2), "gaussian"),
    control.family = list(
        list(),
        list(),
        list(hyper = list(prec = list(initial = 20, fixed = TRUE)))),
    data = data3,
    control.mode = list(
        theta = c(3, 1.6, 0),
        restart = TRUE),
    control.predictor = list(link = lnk),
    verbose = TRUE##, inla.call = 'remote'
)
