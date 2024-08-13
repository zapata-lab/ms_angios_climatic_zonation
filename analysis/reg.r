
library(data.table)
library(INLA)
library(ape)
library(MCMCglmm)

# read data
d = readRDS("/Users/quintero/repos/ms_angios_climatic_zonation/final_data_table.rds")

# set as data.table
setDT(d)

# read tree
tree = read.tree('/Users/quintero/repos/ms_angios_climatic_zonation/tree.tre')
tree$node.label = NULL

# precision matrix
ivcv = MCMCglmm::inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv

# match with data
d[,sp1 := gsub(' ', '_', sp1_species)]
d[,sp2 := gsub(' ', '_', sp2_species)]

d[,spid1 := match(d[,sp1], rownames(ivcv))]
d[,spid2 := match(d[,sp2], rownames(ivcv))]


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# explain range
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# organize for each species
d[, pair := 1:nrow(d)]

# 1st
d1 = d[, .(sp1_species, sp1_n_collections, sp1_mean_latitude, sp1_elev_range,
           sp1_temp_range, sp1_precip_range, region, pair, spid1)]

setnames(d1, c('sp', 'n', 'l', 'er', 'tr', 'pr', 'r', 'pid', 'spid'))

# 2st
d2 = d[, .(sp2_species, sp2_n_collections, sp2_mean_latitude, sp2_elev_range,
           sp2_temp_range, sp2_precip_range, region, pair, spid2)]
setnames(d2, c('sp', 'n', 'l', 'er', 'tr', 'pr', 'r', 'pid', 'spid'))

di = rbind(d1,d2)


# define random effects
di[, `:=`(b_0 = pid,
          b_1 = pid,
          b_2 = pid)]
np = nrow(d)

# make combined tropical and extra-tropical
di[r == 'Tropical', rc := "T"]
di[r == 'N. Temperate' | r == 'S. Temperate', rc := "E"]

# scale variables
di[, n_s  := scale(n)]
di[, l_s  := scale(l)]


# complexity prior
pcprior = list(prec = list(
  prior="pc.prec",
  param = c(20, 0.1)) 
)

#########
# elevation regression on the tropics/extra-tropics level
#########
f = er ~ n_s + rc + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
r = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(r)

# make posterior marginal effect plot
eTe = r$marginals.fixed$`(Intercept)`
eTr = r$marginals.fixed$rc
eTr[,1] = eTe[,1] + eTr[,1]

x = seq(min(eTe[,1]), max(eTe[,1]), 0.001)
y = inla.dmarginal(x, eTe)
eTe = matrix(c(x,y), nrow = length(x))
x = seq(min(eTr[,1]), max(eTr[,1]), 0.001)
y = inla.dmarginal(x, eTr)
eTr = matrix(c(x,y), nrow = length(x))

# temperature regression on the tropics/extra-tropics level
f = tr ~ n_s + rc + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
r = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(r)

# make posterior marginal effect plot
tTe = r$marginals.fixed$`(Intercept)`
tTr = r$marginals.fixed$rc
tTr[,1] = tTe[,1] + tTr[,1]

x = seq(min(tTe[,1]), max(tTe[,1]), 0.001)
y = inla.dmarginal(x, tTe)
tTe = matrix(c(x,y), nrow = length(x))
x = seq(min(tTr[,1]), max(tTr[,1]), 0.001)
y = inla.dmarginal(x, tTr)
tTr = matrix(c(x,y), nrow = length(x))


pdf('~/data/ms_angios_climatic_zonation/plots/range_trop_temp.pdf', 
  height = 6, width = 8)
  par(mfrow = c(1,2))

  plot(1, type = 'n', xlim = c(-0.5,1.5), ylim = c(range(di[,er])), 
    xlab = '', ylab = 'Elevational range', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:1, labels = c('Tropical', 'Temperate'))

  y = di[rc == 'T', er]
  points((runif(length(y)) - 0.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[rc == 'E', er]
  points((runif(length(y)) + 1.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*eTr[,2]/max(eTr[,2]),exp(eTr[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*eTe[,2]/max(eTe[,2]) + 1,exp(eTe[,1]), col = '#02315EEE', border = NA)


  plot(1, type = 'n', xlim = c(-0.5,1.5), ylim = c(range(di[,tr])), 
    xlab = '', ylab = 'Temperature range', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:1, labels = c('Tropical', 'Temperate'))

  y = di[rc == 'T', tr]
  points((runif(length(y)) - 0.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[rc == 'E', tr]
  points((runif(length(y)) + 1.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*tTr[,2]/max(tTr[,2]),exp(tTr[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*tTe[,2]/max(tTe[,2]) + 1,exp(tTe[,1]), col = '#02315EEE', border = NA)

dev.off()


#########
# elevation regression on the S/T/N latitudes level
#########
f = er ~ n_s + r + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
r = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(r)

# make posterior marginal effect plot
eN = r$marginals.fixed$`(Intercept)`
eS = r$marginals.fixed$`rS. Temperate`
eT = r$marginals.fixed$rTropical
eS[,1] = eN[,1] + eS[,1]
eT[,1] = eN[,1] + eT[,1]


x = seq(min(eN[,1]), max(eN[,1]), 0.001)
y = inla.dmarginal(x, eN)
eN = matrix(c(x,y), nrow = length(x))
x = seq(min(eS[,1]), max(eS[,1]), 0.001)
y = inla.dmarginal(x, eS)
eS = matrix(c(x,y), nrow = length(x))
x = seq(min(eT[,1]), max(eT[,1]), 0.001)
y = inla.dmarginal(x, eT)
eT = matrix(c(x,y), nrow = length(x))


# temperature regression on the S/T/N latitudes level
f = tr ~ n_s + r + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
r = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(r)

# make posterior marginal effect plot
tN = r$marginals.fixed$`(Intercept)`
tS = r$marginals.fixed$`rS. Temperate`
tT = r$marginals.fixed$rTropical
tS[,1] = tN[,1] + tS[,1]
tT[,1] = tN[,1] + tT[,1]


x = seq(min(tN[,1]), max(tN[,1]), 0.001)
y = inla.dmarginal(x, tN)
tN = matrix(c(x,y), nrow = length(x))
x = seq(min(tS[,1]), max(tS[,1]), 0.001)
y = inla.dmarginal(x, tS)
tS = matrix(c(x,y), nrow = length(x))
x = seq(min(tT[,1]), max(tT[,1]), 0.001)
y = inla.dmarginal(x, tT)
tT = matrix(c(x,y), nrow = length(x))


pdf('~/data/ms_angios_climatic_zonation/plots/range_S_T_N.pdf', 
  height = 6, width = 8)
  par(mfrow = c(1,2))

  plot(1, type = 'n', xlim = c(-0.5,2.5), ylim = c(range(di[,er])), 
    xlab = '', ylab = 'Elevational range', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:2, labels = c('Southern','Tropical', 'Northern'))

  y = di[r == 'S. Temperate', er]
  points(runif(length(y))*0.5 - 0.25, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[r == 'Tropical', er]
  points(runif(length(y))*0.5 + 0.75, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[r == 'N. Temperate', er]
  points(runif(length(y))*0.5 + 1.75, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*eS[,2]/max(eS[,2]),     exp(eS[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*eT[,2]/max(eT[,2]) + 1, exp(eT[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*eN[,2]/max(eN[,2]) + 2, exp(eN[,1]), col = '#02315EEE', border = NA)


  plot(1, type = 'n', xlim = c(-0.5,2.5), ylim = c(range(di[,tr])), 
    xlab = '', ylab = 'Temperature range', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:2, labels = c('Southern','Tropical', 'Northern'))

  y = di[r == 'S. Temperate', tr]
  points(runif(length(y))*0.5 - 0.25, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[r == 'Tropical', tr]
  points(runif(length(y))*0.5 + 0.75, y, pch = 20, col = rgb(0,0,0,0.03))
  y = di[r == 'N. Temperate', tr]
  points(runif(length(y))*0.5 + 1.75, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*tS[,2]/max(tS[,2]),     exp(tS[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*tT[,2]/max(tT[,2]) + 1, exp(tT[,1]), col = '#02315EEE', border = NA)
  polygon(0.3*tN[,2]/max(tN[,2]) + 2, exp(tN[,1]), col = '#02315EEE', border = NA)

dev.off()


#########
# regression with latitude
#########
f = er ~ n_s + l_s + I(l_s^2) + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
re = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(re)

f = tr ~ n_s + l_s + I(l_s^2) + f(spid, model = "generic0", Cmatrix = ivcv, hyper = pcprior)
rt = inla(f, 
  data = di, family = "lognormal", 
  control.compute=list(return.marginals.predictor=TRUE))
summary(rt)


pdf('~/data/ms_angios_climatic_zonation/plots/range_lat.pdf', 
  height = 6, width = 8)

  par(mfrow = c(1,2))

  plot(di[,l_s], di[,er], xlab = 'Latitude', 
    ylab = 'Elevational range', xaxt = 'n', las =  1, bty = 'n', pch = 20,
    col = rgb(0,0,0,0.03))
  x = seq(min(di[,l_s]), max(di[,l_s]), 0.01)
  c = attr(scale(di[,l]), "scaled:center")
  s = attr(scale(di[,l]), "scaled:scale")
  rl = round(range(di[,l]))
  labs = seq(rl[1], rl[2], 15)
  axis(1, at = (labs - c)/s, labels = labs)

  b_0 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$`(Intercept)`)
  b_1 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$l_s)
  b_2 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$`I(l_s^2)`)
  y = exp(mean(b_0[1,]) +  mean(b_1[1,])*x + mean(b_2[1,])*x^2)
  lines(x, y, col = '#02315E', lwd = 2)
  y = exp(mean(b_0[2,1]) +  mean(b_1[2,1])*x + mean(b_2[2,1])*x^2)
  lines(x, y, col = '#02315E', lty = 2)
  y = exp(mean(b_0[2,2]) +  mean(b_1[2,2])*x + mean(b_2[2,2])*x^2)
  lines(x, y, col = '#02315E', lty = 2)

  plot(di[,l_s], di[,tr], xlab = 'Latitude', 
    ylab = 'Temperature range', xaxt = 'n', las =  1, bty = 'n', pch = 20,
    col = rgb(0,0,0,0.03))
  x = seq(min(di[,l_s]), max(di[,l_s]), 0.01)
  c = attr(scale(di[,l]), "scaled:center")
  s = attr(scale(di[,l]), "scaled:scale")
  rl = round(range(di[,l]))
  labs = seq(rl[1], rl[2], 15)
  axis(1, at = (labs - c)/s, labels = labs)

  b_0 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$`(Intercept)`)
  b_1 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$l_s)
  b_2 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$`I(l_s^2)`)
  y = exp(mean(b_0[1,]) +  mean(b_1[1,])*x + mean(b_2[1,])*x^2)
  lines(x, y, col = '#02315E', lwd = 2)
  y = exp(mean(b_0[2,1]) +  mean(b_1[2,1])*x + mean(b_2[2,1])*x^2)
  lines(x, y, col = '#02315E', lty = 2)
  y = exp(mean(b_0[2,2]) +  mean(b_1[2,2])*x + mean(b_2[2,2])*x^2)
  lines(x, y, col = '#02315E', lty = 2)
dev.off()





#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# explain overlap
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# scale temperature overlap
d[,eo := elevation_overlap]
d[,to := temperature_overlap/12.0]

# combine regions
d[region == 'Tropical', rc := "T"]
d[region == 'N. Temperate' | region == 'S. Temperate', rc := "E"]

# min number of occurrences
d[,n_mn := scale(pmin(sp1_n_collections, sp2_n_collections))]

# rearrange tree to only contain the pendant edge to the sister species
t0 = keep.tip(tree, d[,sp1])

# remove branch lengths of sister species from the tree
t0$tip.label
d = d[match(t0$tip.label,sp1)]
el = d[, pair_age]

t0$edge.length[t0$edge[,2] %in% 1:743] = 
  t0$edge.length[t0$edge[,2] %in% 1:743] - el

# scale tree
s = branching.times(tree)[1] - min(el)

t0$edge.length = t0$edge.length/s

# precision matrix
ivcv0 = MCMCglmm::inverseA(t0, nodes = "TIPS", scale = F)$Ainv

d[,spid := match(d[,sp1], rownames(ivcv0))]



#########
# regression on the tropics/extra-tropics level
#########

# elevation overlap
f = eo ~ pair_age + n_mn + rc + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior)
re = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(re)

# make posterior marginal effect plot
eTe = re$marginals.fixed$`(Intercept)`
eTr = re$marginals.fixed$rcT
eTr[,1] = eTe[,1] + eTr[,1]

x = seq(min(eTe[,1]), max(eTe[,1]), 0.001)
y = inla.dmarginal(x, eTe)
eTe = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))

x = seq(min(eTr[,1]), max(eTr[,1]), 0.001)
y = inla.dmarginal(x, eTr)
eTr = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))


# temperature overlap
f = to ~ pair_age + n_mn + rc + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior) 
rt = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(rt)

# make posterior marginal effect plot
tTe = rt$marginals.fixed$`(Intercept)`
tTr = rt$marginals.fixed$rcT
tTr[,1] = tTe[,1] + tTr[,1]

x = seq(min(tTe[,1]), max(tTe[,1]), 0.001)
y = inla.dmarginal(x, tTe)
tTe = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))

x = seq(min(tTr[,1]), max(tTr[,1]), 0.001)
y = inla.dmarginal(x, tTr)
tTr = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))

pdf('~/data/ms_angios_climatic_zonation/plots/over_trop_temp.pdf', 
  height = 6, width = 8)
  par(mfrow = c(1,2))

  plot(1, type = 'n', xlim = c(-0.5,1.5), ylim = c(range(d[,eo])), 
    xlab = '', ylab = 'Elevational overlap', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:1, labels = c('Tropical', 'Temperate'))

  y = d[rc == 'T', eo]
  points((runif(length(y)) - 0.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[rc == 'E', eo]
  points((runif(length(y)) + 1.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*eTr[,2]/max(eTr[,2]),    eTr[,1], col = '#02315EEE', border = NA)
  polygon(0.3*eTe[,2]/max(eTe[,2]) + 1,eTe[,1], col = '#02315EEE', border = NA)


  plot(1, type = 'n', xlim = c(-0.5,1.5), ylim = c(range(d[,to])), 
    xlab = '', ylab = 'Temperature overlap', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:1, labels = c('Tropical', 'Temperate'))

  y = d[rc == 'T', to]
  points((runif(length(y)) - 0.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[rc == 'E', to]
  points((runif(length(y)) + 1.5)*0.5, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*tTr[,2]/max(tTr[,2]),    tTr[,1], col = '#02315EEE', border = NA)
  polygon(0.3*tTe[,2]/max(tTe[,2]) + 1,tTe[,1], col = '#02315EEE', border = NA)
dev.off()





#########
# regression on the S, T, N level
#########

# elevation overlap
f = eo ~ pair_age + n_mn + region + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior) 
re = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(re)

# make posterior marginal effect plot
eN = re$marginals.fixed$`(Intercept)`
eT = re$marginals.fixed$regionTropical
eS = re$marginals.fixed$`regionS. Temperate`
eT[,1] = eT[,1] + eN[,1]
eS[,1] = eS[,1] + eN[,1]

x = seq(min(eT[,1]), max(eT[,1]), 0.001)
y = inla.dmarginal(x, eT)
eT = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))
x = seq(min(eS[,1]), max(eS[,1]), 0.001)
y = inla.dmarginal(x, eS)
eS = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))
x = seq(min(eN[,1]), max(eN[,1]), 0.001)
y = inla.dmarginal(x, eN)
eN = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))


# temperature overlap
f = to ~ pair_age + n_mn + region + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior) 
rt = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(rt)

# make posterior marginal effect plot
tN = rt$marginals.fixed$`(Intercept)`
tT = rt$marginals.fixed$regionTropical
tS = rt$marginals.fixed$`regionS. Temperate`
tT[,1] = tT[,1] + tN[,1]
tS[,1] = tS[,1] + tN[,1]

x = seq(min(tT[,1]), max(tT[,1]), 0.001)
y = inla.dmarginal(x, tT)
tT = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))
x = seq(min(tS[,1]), max(tS[,1]), 0.001)
y = inla.dmarginal(x, tS)
tS = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))
x = seq(min(tN[,1]), max(tN[,1]), 0.001)
y = inla.dmarginal(x, tN)
tN = matrix(c(exp(x)/(1 + exp(x)),y), nrow = length(x))

pdf('~/data/ms_angios_climatic_zonation/plots/over_S_T_N.pdf', 
  height = 6, width = 8)
  par(mfrow = c(1,2))

  plot(1, type = 'n', xlim = c(-0.5,2.5), ylim = c(range(d[,eo])), 
    xlab = '', ylab = 'Elevational overlap', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:2, labels = c('Southern','Tropical', 'Northern'))

  y = d[region == 'S. Temperate', eo]
  points(runif(length(y))*0.5 - 0.25, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[region == 'Tropical', eo]
  points(runif(length(y))*0.5 + 0.75, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[region == 'N. Temperate', eo]
  points(runif(length(y))*0.5 + 1.75, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*eS[,2]/max(eS[,2]),     eS[,1], col = '#02315EEE', border = NA)
  polygon(0.3*eT[,2]/max(eT[,2]) + 1, eT[,1], col = '#02315EEE', border = NA)
  polygon(0.3*eN[,2]/max(eN[,2]) + 2, eN[,1], col = '#02315EEE', border = NA)


  plot(1, type = 'n', xlim = c(-0.5,2.5), ylim = c(range(d[,to])), 
    xlab = '', ylab = 'Temperature overlap', xaxt = 'n', las =  1, bty = 'n')
  axis(1, at = 0:2, labels = c('Southern','Tropical', 'Northern'))

  y = d[region == 'S. Temperate', to]
  points(runif(length(y))*0.5 - 0.25, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[region == 'Tropical', to]
  points(runif(length(y))*0.5 + 0.75, y, pch = 20, col = rgb(0,0,0,0.03))
  y = d[region == 'N. Temperate', to]
  points(runif(length(y))*0.5 + 1.75, y, pch = 20, col = rgb(0,0,0,0.03))

  polygon(0.3*tS[,2]/max(tS[,2]),     tS[,1], col = '#02315EEE', border = NA)
  polygon(0.3*tT[,2]/max(tT[,2]) + 1, tT[,1], col = '#02315EEE', border = NA)
  polygon(0.3*tN[,2]/max(tN[,2]) + 2, tN[,1], col = '#02315EEE', border = NA)

dev.off()





#########
# regression with latitude
#########
d[, l   := (sp1_mean_latitude + sp2_mean_latitude)/2]
d[, l_s := scale(l)]

# elevation overlap
f = eo ~ pair_age + n_mn + l_s + I(l_s^2) + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior) 
re = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(re)

# temperature overlap
f = to ~ pair_age + n_mn + l_s + I(l_s^2) + f(spid, model = "generic0", Cmatrix = ivcv0, hyper = pcprior) 
rt = inla(f, data = d, family = "beta", 
  control.family  = list(beta.censor.value = 0.00001),
  control.compute = list(return.marginals.predictor=TRUE))
summary(rt)




pdf('~/data/ms_angios_climatic_zonation/plots/over_lat.pdf', 
  height = 6, width = 8)

  par(mfrow = c(1,2))

  plot(d[,l_s], d[,eo], xlab = 'Latitude', 
    ylab = 'Elevational overlap', xaxt = 'n', las =  1, bty = 'n', pch = 20,
    col = rgb(0,0,0,0.03))
  x = seq(min(d[,l_s]), max(d[,l_s]), 0.01)
  c = attr(scale(d[,l]), "scaled:center")
  s = attr(scale(d[,l]), "scaled:scale")
  rl = round(range(d[,l]))
  labs = seq(rl[1], rl[2], 15)
  axis(1, at = (labs - c)/s, labels = labs)

  b_0 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$`(Intercept)`)
  b_1 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$l_s)
  b_2 = inla.hpdmarginal(c(0.025, 0.975), re$marginals.fixed$`I(l_s^2)`)
  y = mean(b_0[1,]) +  mean(b_1[1,])*x + mean(b_2[1,])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lwd = 2)
  y = mean(b_0[2,1]) +  mean(b_1[2,1])*x + mean(b_2[2,1])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lty = 2)
  y = mean(b_0[2,2]) +  mean(b_1[2,2])*x + mean(b_2[2,2])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lty = 2)

  plot(d[,l_s], d[,to], xlab = 'Latitude', 
    ylab = 'Temperature overlap', xaxt = 'n', las =  1, bty = 'n', pch = 20,
    col = rgb(0,0,0,0.03))
  x = seq(min(d[,l_s]), max(d[,l_s]), 0.01)
  c = attr(scale(d[,l]), "scaled:center")
  s = attr(scale(d[,l]), "scaled:scale")
  rl = round(range(d[,l]))
  labs = seq(rl[1], rl[2], 15)
  axis(1, at = (labs - c)/s, labels = labs)

  b_0 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$`(Intercept)`)
  b_1 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$l_s)
  b_2 = inla.hpdmarginal(c(0.025, 0.975), rt$marginals.fixed$`I(l_s^2)`)
  y = mean(b_0[1,]) +  mean(b_1[1,])*x + mean(b_2[1,])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lwd = 2)
  y = mean(b_0[2,1]) +  mean(b_1[2,1])*x + mean(b_2[2,1])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lty = 2)
  y = mean(b_0[2,2]) +  mean(b_1[2,2])*x + mean(b_2[2,2])*x^2
  lines(x, exp(y)/(1+exp(y)), col = '#02315E', lty = 2)
dev.off()




library(zoib)

d[, al_s := scale(abs(l))]

# run analyses
r_zoib_e = zoib(eo ~ pair_age + n_mn + al_s | pair_age + n_mn + al_s | al_s | al_s, 
  data = d, random = 0, zero.inflation = TRUE, one.inflation = TRUE, 
  n.iter = 1050, n.thin = 5, n.burn=50)
summary(r_zoib_e$coeff)

# make new data to predict effect of latitude while holding pari age and number of records
#constant
xnew  = data.table(pair_age = mean(d[,pair_age]),
                   al_s = seq(-1, 1, 0.1), 
                   n_mn = mean(d[,n_mn]))

r_zoib_t = zoib(to ~ pair_age + n_mn + al_s | pair_age + n_mn + al_s | al_s, 
  data = d, random = 0, zero.inflation = TRUE, one.inflation = FALSE, 
  n.iter = 1050, n.thin = 5, n.burn=50)
summary(r_zoib_t$coeff)


xls = range(d[,abs(l)])
xls = signif(seq(xls[1], xls[2], length.out = 5), 2)


# make prediction plots
pdf('~/data/ms_angios_climatic_zonation/plots/over_lat_zoib.pdf', 
  height = 5, width = 8)

  par(mfrow = c(1,2), bty = 'n', las = 1)

  # predict effect of latitude
  ypred = pred.zoib(r_zoib_e, xnew)
  ynew = ypred$summary

  plot(xnew$al_s, ynew[,2], xlab = 'Absolute latitude', type = 'l', 
    ylim = range(ynew[,7:8]), ylab = 'Posterior predicted elevation overlap', 
    lwd = 2, xaxt = 'n', col = '#02315E')
  axis(1, at = seq(-1, 1, length.out = 5), labels = xls)
  lines(xnew$al_s, ynew[,7], lty = 2, lwd = 2, col = '#02315E')
  lines(xnew$al_s, ynew[,8], lty = 2, lwd = 2, col = '#02315E')

  # predict effect of latitude
  ypred = pred.zoib(r_zoib_t, xnew)
  ynew = ypred$summary

  plot(xnew$al_s, ynew[,2], xlab = 'Absolute latitude', type = 'l', 
    ylim = range(ynew[,7:8]), ylab = 'Posterior predicted temperature overlap',
    lwd = 2, xaxt = 'n', col = '#02315E')
  axis(1, at = seq(-1, 1, length.out = 5), labels = xls)
  lines(xnew$al_s, ynew[,7], lty = 2, lwd = 2, col = '#02315E')
  lines(xnew$al_s, ynew[,8], lty = 2, lwd = 2, col = '#02315E')


dev.off()


