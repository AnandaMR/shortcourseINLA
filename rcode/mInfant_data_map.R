library(ggplot2)
    

setwd(here::here(""))
getwd()

nasc <- read.csv2(
    file = "data/nasc.csv", skip = 3,
    nrows = 855, encoding = 'latin1',
    na.strings = '-')

head(nasc,2)
tail(nasc,2)

iobt <- read.csv2(
    file = "data/iobt.csv", skip = 3,
    nrows = 855, encoding = 'latin1',
    na.strings = '-')

head(iobt,2)
tail(iobt,2)

(jj <- 2:(ncol(iobt)-1))
years <- as.integer(substr(
    colnames(iobt)[jj], 2, 5))
years

if(!file.exists("data/mInfant.rds")) {
    muns <- c("BETIM", "BELO HORIZONTE", "CONTAGEM",
              "NOVA LIMA", "RIBEIRAO DAS NEVES", "SABARA")
    (nm <- length(muns))
    (i1 <- sapply(muns, grep, nasc$Mun))
    (i2 <- sapply(muns, grep, iobt$Mun))
    (nyears <- length(years))
    
    mInfant <- data.frame(
        year = rep(years, nm),
        mun = rep(muns, each = nyears),
        nasc = as.vector(t(nasc[, jj])[, i1]),
        iobt = as.vector(t(iobt[, jj])[, i2])
    )
    saveRDS(mInfant, "data/mInfant.rds")
    
    str(mInfant)
    
    ggplot(mInfant) + theme_minimal() +
        geom_line(aes(x = year, y = iobt/nasc, color = mun))

}

## package to work work with spatial/maps

library(sf)

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

head(mg_mun)

o1 <- pmatch(substr(nasc$Mun, 1, 6),
             substr(mg_mun$code, 1, 6))
summary(o1)
nasc[is.na(o1),]

o2 <- pmatch(substr(iobt$Mun, 1, 6),
             substr(mg_mun$code, 1, 6))
summary(o2)
nasc[is.na(o2),]

names(nasc)

mg_mun$nasc <- nasc[o1[complete.cases(o1)], 30]
mg_mun$iobt <- iobt[o2[complete.cases(o2)], 30]

head(mg_mun, 3)

summary(mg_mun$nasc)

bb <- st_bbox(mg_mun)
bb

rr <- c(x = diff(bb[c(1,3)]),
        y = diff(bb[c(2,4)]))
rr

rr * 200

png("iobt_map_mg.png", 2300, 1800, res = 300)
ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = iobt/nasc), color = 'transparent') +
    scale_fill_distiller(
        palette = "RdBu", direction = -1
    )
dev.off()

system("eog iobt_map_mg.png &")


library(spdep)

viz <- poly2nb(mg_mun)

viz

xy <- st_coordinates(st_centroid(mg_mun))

par(mar = c(0,0,0,0))
plot(viz, xy)

nb2INLA("grafo", viz)

## modelo y_i ~ Poisson(r_i * N_i)

ff1 <- iobt ~
    f(area, model = "besag", graph = "grafo",
      scale.model = TRUE)

library(INLA)

mg_mun$area <- 1:nrow(mg_mun)

fit1 <- inla(
    formula = ff1,
    family = "poisson",
    data = st_drop_geometry(mg_mun)[c("nasc", "iobt", "area")],
    E = nasc,
    control.predictor = list(link = 1)
)

fit1$summary.fixed

exp(-4.09)

plot(fit1$marginals.random$area[[1]], type = "l", lwd = 5, xlim = c(-1,1))
lines(fit1$marginals.random$area[[2]], col = 2, lwd = 5)
lines(fit1$marginals.random$area[[3]], col = 3, lwd = 5)
lines(fit1$marginals.random$area[[4]], col = 4, lwd = 5)

mg_mun$log_ri <- fit1$summary.random$area$mean

ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = log_ri)) +
    scale_fill_distiller(
        palette = "RdBu"
    )

### lambda0 = sum_i obitos_i / sum_i nascidos_i
lambda0 <- sum(mg_mun$iobt, na.rm = TRUE) / sum(mg_mun$nasc)
lambda0

mg_mun$E <- lambda0 * mg_mun$nasc

### \lambda = E * exp( preditor linear )

fit2 <- inla(
    formula = ff1,
    family = "poisson",
    data = st_drop_geometry(mg_mun)[c("nasc", "iobt", "area", "E")],
    E = E,
    control.predictor = list(link = 1)
)

fit2$summary.fixed

mg_mun$log_smr <- fit2$summary.random$area$mean

ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = log_smr)) +
    scale_fill_distiller(
        palette = "RdBu"
    )

names(inla.models())
names(inla.models()$latent)
## inla.models()$latent$besagproper
inla.doc("besagproper")


ff3 <- iobt ~
    f(area, model = "besagproper", graph = "grafo")

fit3 <- inla(
    formula = ff3,
    family = "poisson",
    data = st_drop_geometry(mg_mun)[c("nasc", "iobt", "area", "E")],
    E = E,
    control.predictor = list(link = 1)
)

fit3$summary.fixed

mg_mun$log_smr3 <- fit3$summary.random$area$mean

ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = log_smr3)) +
    scale_fill_distiller(
        palette = "RdBu"
    )


fit2$summary.hyperpar
fit3$summary.hyperpar

## bym2 :  ( sqrt(\phi) * u_i + sqrt(1- \phi) * v_i) / sqrt(tau)
##  \phi : quao forte eh a dependencia espacial
ff4 <- iobt ~
    f(area, model = "bym2", graph = "grafo")

fit4 <- inla(
    formula = ff4,
    family = "poisson",
    data = st_drop_geometry(mg_mun)[c("nasc", "iobt", "area", "E")],
    E = E,
    control.predictor = list(link = 1)
)

fit4$summary.hyper

mg_mun$log_smr4 <- fit3$summary.random$area$mean[1:nrow(mg_mun)]

ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = log_smr4)) +
    scale_fill_distiller(
        palette = "RdBu"
    )

