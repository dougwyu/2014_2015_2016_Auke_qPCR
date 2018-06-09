library(ggplot2)


## Upstream Rating Curve ####

load("upstream.RData")


dat.u <- up[, c(3, 2)]

min.rss <- function(data, par) {
  
  with(data, sum((cfs - (par[3] * ((gauge - par[1])^par[2])))^2))
  
}

load("parupstream.RData")

op.u <- optim(par = par.u, fn = min.rss, data = dat.u)

par.u <- op.u$par

save(par.u, file = "parupstream.RData")

gauge.u.predict <- seq(20.33, 23.01, 0.02)

cfs.u.predict <- par.u[3] * ((gauge.u.predict - par.u[1])^par.u[2])



ggplot() + geom_line(aes(gauge.u.predict, cfs.u.predict)) + geom_point(data = up, 
  aes(gauge, cfs), shape = 1)




## QFS Function ####


cfs <- function(gaugeheight) {
  
  
  par.u[3] * ((gaugeheight - par.u[1])^par.u[2])
  
}

cfs(22)




