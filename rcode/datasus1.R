## remotes::install_github("rfsaldanha/microdatasus")

library(microdatasus)

sinasc_all <- fetch_datasus(
  year_start = 2024, 
  year_end = 2024, 
  uf = "MG", 
  information_system = "SINASC"
)

dim(sinasc_all)

names(sinasc_all)

table(sinasc_all$CONSULTAS)

table(sinasc_all$APGAR1)
table(as.numeric(sinasc_all$APGAR1)<=6)

summary(as.numeric(sinasc_all$PESO)/1000)


tbmuni <- sort(unique(sinasc_all$CODMUNRES))
length(tbmuni)
head(tbmuni)

dataf <- data.frame(
    peso = as.numeric(sinasc_all$PESO)/1000,
    ap1b = (as.numeric(sinasc_all$APGAR1)<=6)+0L,
    mun  = pmatch(sinasc_all$CODMUNRES,
                  tbmuni,
                  duplicates.ok = TRUE)
)

str(dataf)

summary(dataf)

m0_ml <- glm(
    formula = ap1b ~ peso,
    family = binomial,
    data = dataf)

nrow(dataf)
str(m0_ml$fitted)

library(INLA)

m0_b <- inla(
    formula = ap1b ~ peso,
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))

coef(summary(m0_ml))
m0_b$summary.fixed

## m1: eta_i = beta_0 + beta_1 * peso_i + u_mun(i)
##     P(y_i=1) = p_i = 1/(1+exp(-eta_i))
m1_b <- inla(
    formula = ap1b ~ peso + f(mun, model = 'iid'),
    family = "binomial",
    data = dataf,
    control.predictor = list(link = 1))

m1_b$summary.fixed

dim(m1_b$summary.random$mun)

m1_b$summary.hy

1/sqrt(6.9)

p2 <- 1/(1+exp(-(-0.030 -0.93669 * 2)))
p2

p3 <- 1/(1+exp(-(-0.030 -0.93669 * 3)))
p3

p2/p3

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


head(m1_b$summary.random$mun,3)

omuni <- pmatch(tbmuni[-1], substr(mg_muni$code_muni,1,6))
stopifnot(!any(is.na(omuni)))

summary(m1_b$summary.random$mun$mean[-1][omuni])

mg_muni$log_r <- m1_b$summary.random$mun$mean[-1][omuni]

ggplot() +
    geom_sf(data = mg_muni,
            aes(fill = log_r),
            color = "transparent", size = 0.15) +
  labs(title = "Municipalities of Minas Gerais",
       caption = "Source: IBGE via geobr") +
    theme_minimal() +
    scale_fill_distiller(
        palette = "RdBu"
    )

which(m1_b$summary.random$mun$mean>2)

pmatch(tbmuni[704], substr(mg_muni$code_muni,1,6))
mg_muni[703,]

dataf[dataf$mun==704,]
