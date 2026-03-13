## remotes::install_github("rfsaldanha/microdatasus")

library(microdatasus)

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

table(sinasc$CONSULTAS)

table(sinasc$APGAR1)
table(as.numeric(sinasc$APGAR1)<=6)

summary(as.numeric(sinasc$PESO)/1000)

tbmuni <- sort(unique(sinasc$CODMUNRES))
length(tbmuni)
head(tbmuni)

table(sinasc$APGAR1,
      sinasc$APGAR5)

dataf <- data.frame(
    peso = as.numeric(sinasc$PESO)/1000,
    ap1b = (as.numeric(sinasc$APGAR1)<=6)+0L,
    ap5b = (as.numeric(sinasc$APGAR5)<=6)+0L,
    mun  = pmatch(sinasc$CODMUNRES,
                  tbmuni,
                  duplicates.ok = TRUE)
)

str(dataf)

summary(dataf)

library(INLA)

fit0a1 <- inla(
    formula = ap1b ~ peso,
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))
fit0a5 <- inla(
    formula = ap5b ~ peso,
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))

fit0a1$summary.fixed
fit0a5$summary.fixed

## m1: eta_i = beta_0 + beta_1 * peso_i + u_mun(i)
##     P(y_i=1) = p_i = 1/(1+exp(-eta_i))
fit1a1 <- inla(
    formula = ap1b ~ peso + f(mun, model = 'iid'),
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))
fit1a5 <- inla(
    formula = ap5b ~ peso + f(mun, model = 'iid'),
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))

fit1a1$summary.fixed
fit1a5$summary.fixed

fit1a1$summary.hyperpar
fit1a5$summary.hyperpar

n <- nrow(dataf)
dat2 <- data.frame(
    ap1 = c(dataf$ap1, rep(NA, n)),
    ap5 = c(rep(NA, n), dataf$ap5),
    b0_a1 = c(rep(1., n), rep(NA, n)), ##rep(1:0, c(n,n)),
    b0_a5 = c(rep(NA, n), rep(1., n)), ##rep(0:1, c(n,n)),
    peso_a1 = c(dataf$peso, rep(NA, n)),
    peso_a5 = c(rep(NA, n), dataf$peso)
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

dat2$mun1 <- c(dataf$mun, rep(NA,n))
dat2$mun2 <- c(rep(NA,n), dataf$mun)

fit1_b <- inla(
    formula = list(ap1, ap5) ~ -1 +
        b0_a1 + b0_a5 + peso_a1 + peso_a5 +
        f(mun1, model = "iid") + f(mun2, model = "iid"),
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
        f(mun1, model = "iid") + f(mun2, copy = 'mun1', fixed = FALSE),

fit2_b <- inla(
    formula = ff2b
    family = rep("binomial", 2),
    data = dat2,
    control.predictor = list(link = 1)
)

fit2_b$summary.hyperpar

##
library(geobr)

mg_muni <- read_municipality(
  code_muni = "MG", 
  year = 2020, 
  showProgress = TRUE
)

head(mg_muni)

dim(mg_muni)

library(ggplot2)

ggplot() +
    geom_sf(data = mg_muni, fill = "#2D3E50",
            color = "#FEBF57", size = 0.15) +
  labs(title = "Municipalities of Minas Gerais",
       caption = "Source: IBGE via geobr") +
    theme_minimal()

head(sinasc, 2)

str(sinasc$CODMUNRES)

id_map <- pmatch(
    x = sinasc$CODMUNRES,
    table = substr(mg_muni$code_muni, 1,6),
    duplicates.ok = TRUE)
summary(id_map)

dat2$area1 <- c(id_map, rep(NA, n))
dat2$area2 <- c(rep(NA, n), id_map)

library(spdep)

viz <- poly2nb(mg_muni)
nb2INLA(file = "grafo", nb = viz) ## 

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
    ap1 = c(dataf$ap1, rep(NA, n),  rep(NA, n)),
    ap5 = c(rep(NA, n), dataf$ap5,  rep(NA, n)),
    zero= c(rep(NA, n), rep(NA,n),  rep(0, n)),
    lnk = rep(1:3, each = n),
    b0_a1 = rep(c(1,0,0), c(n,n,n)),
    b0_a5 = rep(c(0,1,0), c(n,n,n)),
    peso1 = c(dataf$peso, rep(NA, n), dataf$peso),
    peso5 = c(rep(NA, n), dataf$peso, rep(NA, n)), ## eff. differencial de peso em ap5
    area  = c(id_map,     rep(NA, n), id_map),
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
