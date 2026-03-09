## https://github.com/eliaskrainski/shortcourseINLA

## ----setup
options(digits=4)

library(MASS)
library(INLA)
library(ggplot2)
library(gridExtra)
library(viridisLite)

## linear predictor, fig in the slides
par(mfrow=c(1,2), mar=c(3,3,1,1),
    mgp=c(2,1,0), xaxs='i', yaxs='i')
F1 <- 0:10
plot(F1, 5 + 0.5*F1, type='l',
     xlab='F', ylab='5 + 0,5*F')
F2 <- 0:10
s <- outer(F1, F2, function(a,b)
    5 + 0.6*a -0.3*b)
persp(F1, F2, s, xlab='F1', ylab='F2',
      zlab='5 + 0.6*F1 - 0.3*F2', theta=5)


## dataset 
data("Orange")
head(Orange, 3)
tail(Orange, 3)


## plot 
g.o1 <- ggplot(Orange, aes(y=circumference)) + 
    geom_point(pch=19,
               aes(x=age, group=Tree,
                   col=Tree)) + 
  xlab("Age (days)")+ 
  theme_grey(base_size = 20)
g.o1


## model 1
m1 <- inla(circumference ~ age, data=Orange, 
           control.compute = list(cpo = TRUE, waic = TRUE),
           verbose = !TRUE)

m1$cpu.used

plot(m1, plot.opt.trace = TRUE)

m1$summary.fixed

names(m1)

plot(m1$marginals.fixed[[1]])
plot(m1$marginals.fixed[[2]])

m1$summary.hyperpar[1,,drop=T]

plot(m1$marginals.hyperpar[[1]])

## theta = log(tau)

## marginals 
grid.arrange(
    ggplot(as.data.frame(
        m1$marginals.fixed[[1]])) +
    geom_line(aes(x=x,y=y)) +
    xlab(expression(beta[0])) +
    ylab('Density'), 
    ggplot(as.data.frame(
        m1$marginals.fixed[[2]])) +
    geom_line(aes(x=x,y=y)) +
    xlab(expression(beta[1])) +
    ylab('Density'), 
    ggplot(as.data.frame(
        inla.tmarginal(function(x)
            exp(-x/2), 
            m1$internal.marginals.hyperpar[[1]]))) +
    geom_line(aes(x=x,y=y)) +
    xlab(expression(sigma)) +
    ylab('Density'))


## fitted
b.m1 <- m1$summary.fixed$mean
g.o1 + geom_abline(intercept = b.m1[1], slope = b.m1[2])


## cpo and pit 
Orange$logCPO <- log(m1$cpo$cpo)
Orange$PIT <- m1$cpo$pit 

CPO <- ggplot(Orange) + 
          geom_point(aes(x=age, y=logCPO))+ 
          theme_grey(base_size = 24) + 
          xlab(" ")

PIT <- ggplot(Orange) + 
          geom_point(aes(x=age, y=PIT))+ 
          theme_grey(base_size = 24)

grid.arrange(CPO, PIT, nrow=2)


## model 2
f2 <- circumference ~ 1 +
    f(Tree, age, model='iid') 

m2 <- inla(f2, data=Orange,
           control.compute=list(cpo=TRUE, waic = TRUE),
##           verbose = TRUE,
          # control.mode = list(
           #    theta = c(-4.6, 4),
            #   restart = TRUE
#           )
)
grep("fn", m2$logfile, value = TRUE)

m2$misc$nfunc

plot(m2, F, F, F, F, F, F, F, F, plot.opt.trace = TRUE)

m2$summary.fixed

m2$summary.random$Tree

m1$summary.hyperpar
m2$summary.hyperpar

s2e_1 <- inla.tmarginal(function(x) exp(-x/2), 
                        m1$internal.marginals.hyperpar[[1]])
s2e_2 <- inla.tmarginal(function(x) exp(-x/2), 
                        m2$internal.marginals.hyperpar[[1]])

inla.emarginal(function(x) exp(-x/2),
               m1$internal.marginals.hyperpar[[1]])
inla.emarginal(function(x) exp(-x/2),
               m2$internal.marginals.hyperpar[[1]])

inla.emarginal(function(x) exp(-x/2),
               m2$internal.marginals.hyperpar[[2]])


c(m1=-sum(log(m1$cpo$cpo)),
  m2=-sum(log(m2$cpo$cpo)))

summary(m1)
summary(m2)

c(m1 = m1$waic$waic, m2 = m2$waic$waic)

head(Orange)
m1_ml <- lm(circumference ~ age, Orange)
AIC(m1_ml)

## cpo and pit 
Orange$logCPO2 <- log(m2$cpo$cpo)
Orange$PIT2 <- m2$cpo$pit 
grid.arrange(ggplot(Orange) +
             geom_point(aes(x=age, y=logCPO2)), 
             ggplot(Orange) +
             geom_point(aes(x=age, y=PIT2)))


f3 <- circumference ~ 1 + f(age, model = "rw2")

m3 <- inla(f3, data = Orange,
           control.compute=list(cpo=TRUE, waic = TRUE))


f4 <- circumference ~ 1 + f(age, model = "rw2", replicate = Tree)

m4 <- inla(f4, data = Orange,
           control.compute=list(cpo=TRUE, waic = TRUE))

c(m1=-sum(log(m1$cpo$cpo)),
  m2=-sum(log(m2$cpo$cpo)),
  m3=-sum(log(m3$cpo$cpo)),
  m4=-sum(log(m4$cpo$cpo)))

plot(m4)

