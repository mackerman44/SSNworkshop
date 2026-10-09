library(classInt)
library(viridis)
library(spmodel)

# ------------------------------------------------------------------------------
# What is a statistic?
# ------------------------------------------------------------------------------

# First supose these are some real data
n = 10
y = c(5.54, 2.80, 6.29, 2.85, 1.85, 4.43, 2.26, 5.54, 1.45, 0.98)
# compute it's mean
mean(y)
# a "statistic" is any function of data
sd(y)

# ------------------------------------------------------------------------------
# But "Statistics" is more than just computing statistics
# It is about "Estimation"
# ------------------------------------------------------------------------------
mean(y)
# is an estimator of the mean of the process that generated the data
# or a MODEL
set.seed(101)
n = 100
y = rnorm(n, mean = 10, sd = 2)
mean(y)
y = rnorm(n, mean = 20, sd = 2)
mean(y)
sd(y)
# Estimation leads to powerful uses of even very simple models
par(mar = c(5,5,1,1))
plot((-100:100)/20 + mean(y),
	dnorm((-100:100)/20 + mean(y), mean = mean(y), sd = sd(y)),
	type = 'l', lwd = 2, xlab = 'Value', ylab = 'Probability Distribution',
	cex.lab = 2, cex.axis = 1.5)
# probability of normal distribution
pnorm(25, mean = mean(y), sd = sd(y))
# quantile of a normal distribution
qnorm(0.999, mean = mean(y), sd = sd(y))


# ------------------------------------------------------------------------------
# Even more than Estimation, Statistics is about quantifying uncertainty
# about estimators
# ------------------------------------------------------------------------------

# A confidence interval
set.seed(1001)
n = 1000
y = rnorm(n, mean = 20, sd = 4)
mean(y)
sd(y)
c(mean(y) - 2*sd(y)/sqrt(n), mean(y) + 2*sd(y)/sqrt(n))

# ------------------------------------------------------------------------------
# But what is being estimated?  We need to assume some sort of
# model for the process that generated the data
# ------------------------------------------------------------------------------

set.seed(1007)

x = 1:100
# exactly the same models
y1 = 2 + 0.5*x + rnorm(100, mean = 0, sd = 4)
y2 =  rnorm(100, mean = 2 + 0.5*x, sd = 4)

par(mar = c(5,5,1,1))
plot(x, y1, pch = 19, cex.lab = 2, cex.axis = 1.5)
points(x, y2, pch = 2, cex = 1.5)
lines((1:100), 2 + 0.5*(1:100), lwd = 3)

summary(lm(y1 ~ x))
summary(lm(y2 ~ x))

# ------------------------------------------------------------------------------
# Testing
# Historical Approach
# ------------------------------------------------------------------------------

n = 10
y1 = rnorm(n, mean = 5, sd = 2)
y2 = rnorm(n, mean = 4, sd = 2)
t.test(y1, y2)
d1 = data.frame(y =c(y1, y2), x = 
	as.factor(c(rep(0, times = n), rep(1, times = n))))
str(d1)
summary(lm(y ~ x, data = d1))

# ------------------------------------------------------------------------------
# Power
# ------------------------------------------------------------------------------

n = 10
y1 = rnorm(n, mean = 4.1, sd = 2)
y2 = rnorm(n, mean = 4, sd = 2)
d1 = data.frame(y =c(y1,y2), x = 
	as.factor(c(rep(0, times = n), rep(1, times = n))))
str(d1)
summary(lm(y ~ x, data = d1))

# ------------------------------------------------------------------------------
# Spatial Data -- Distances and Covariance Matrices
# ------------------------------------------------------------------------------

n = 10
set.seed(2003)

#simulate some spatial coordinates
sx = runif(n)
sy = runif(n)
par(mar = c(5,5,1,1))
plot(sx, sy, pch = 19, cex = 2, xlim = c(0,1), ylim = c(0,1),
	cex.lab = 2, cex.axis = 1.5)

# put spatial data into a 2 column matrix
sdat = cbind(sx, sy)
# get distances between all spatial locations as an n x n matrix
dmat = as.matrix(dist(sdat))
dmat
# create a diagonal matrix with just 1's on the diagonal
imat = diag(n)
# create a range of distances from 0 to 1
xrange = (1:1000)/1000
# plot the exponential autocorrelation function
par(mar = c(5,5,1,1))
de = 2 # dependent error, partial sill
range = 0.2 
ie = 0.1 # independent error, nugget effect
plot(xrange, de*exp(-xrange/range), type = 'l', lwd = 2,
	xlim = c(0, 1), ylim = c(0, (ie + de)),
	xlab = 'distance', ylab = 'autocovariance',
	cex.lab = 2, cex.axis = 1.5)
points(0, de + ie, cex = 2, pch = 19)
# create a covariance matrix using distances and exponential autocorrelation
# function
partial_sill = 2 # dependent error (de)
range = 0.2
nugget = 0.1 # independent error (ie)
covmat = partial_sill*exp(-dmat/range) + nugget*imat
covmat

# ------------------------------------------------------------------------------
# Simulate Spatial Data 
# ------------------------------------------------------------------------------

L = chol(covmat)
t(L) %*% L
x = 1:n
y1 = 2 + 0.5*x + t(L) %*% rnorm(n, mean = 0, sd = 1)
d1 = data.frame(y = y1, x = x, sx = sx, sy = sy)

# ------------------------------------------------------------------------------
# Fit Spatial Models 
# ------------------------------------------------------------------------------

summary(splm(y1 ~ x, data = d1, xcoord = sx, ycoord = sy))

# ------------------------------------------------------------------------------
# More Simulations 
# ------------------------------------------------------------------------------

ncr = 40
dat = data.frame(
	sx = as.vector(outer(1:ncr, rep(1, times = ncr)))/ncr,
	sy = as.vector(outer(rep(1, times = ncr), 1:ncr))/ncr)

de = 2 # dependent error, partial sill
range = 0.2 
ie = 0.1 # independent error, nugget effect
# create a range of distances from 0 to 1
xrange = (1:1000)/1000
par(mar = c(5,5,1,1))
plot(xrange, de*exp(-xrange/range), type = 'l', lwd = 2,
	xlim = c(0, 1), ylim = c(0, (ie + de)),
	xlab = 'distance', ylab = 'autocovariance',
	cex.lab = 2, cex.axis = 1.5)
points(0, de + ie, cex = 2, pch = 19)

z = sprnorm(
  spcov_params = spcov_params("exponential", de = de, ie = ie, range = range),
  mean = 0,
  samples = 1,
  data = dat,
  xcoord = sx,
  ycoord = sy)
  
nclass = 12
cip = classIntervals(z, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
par(mar = c(5,5,1,1))
plot(dat$sx, dat$sy, pch = 15, col = cont_colors, cex = 3,
	xlab = 'x', ylab = 'y', cex.lab = 2, cex.axis = 1.5)


# ------------------------------------------------------------------------------
# Prediction 
# ------------------------------------------------------------------------------

n = 100
np = 1000
set.seed(2004)

#simulate some spatial coordinates for observed data
sx = runif(n)
sy = runif(n)
#simulate some spatial coordinates for prediction locations
sxp = runif(np)
syp = runif(np)
par(mar = c(5,5,1,1))
plot(sx, sy, pch = 19, cex = 2, xlim = c(0,1), ylim = c(0,1),
	cex.lab = 2, cex.axis = 1.5)
points(sxp, syp, pch = 1, cex = 1)

# put all spatial data into a 2 column matrix
sdat = cbind(c(sx, sxp), c(sy, syp))

# get distances between all spatial locations as an n x n matrix
dmat = as.matrix(dist(sdat))
dim(dmat)
# create a diagonal matrix with just 1's on the diagonal
imat = diag(n + np)
dim(imat)
# create a covariance matrix using distances and exponential autocorrelation
# function
partial_sill = 2 # dependent error (de)
range = 0.5
nugget = 0.1 # independent error (ie)
covmat = partial_sill*exp(-dmat/range) + nugget*imat
dim(covmat)
# Simulate Spatial Data 
L = chol(covmat)
x = rnorm(n + np)
xo = x[1:n]
xp = x[(n + 1):(n + np)]
y = 2 + 0.5*x + t(L) %*% rnorm(n + np, mean = 0, sd = 1)
yo = y[1:n]
yp = y[(n + 1):(n + np)]
d1 = data.frame(y = yo, x = xo, sx = sx, sy = sy)
# Fit Spatial Models 
splmfit = splm(y ~ x, data = d1, xcoord = sx, ycoord = sy)
summary(splmfit)
# put spatial prediction locations into a data.frame 
p1 = data.frame(sx = sxp, sy = syp, x = xp)
# make the predictions
predfit = predict(splmfit, p1, se.fit = TRUE)
str(predfit)
# plot simulated data versus predicted data
par(mar = c(5,5,1,1))
plot(yp, predfit$fit, pch = 19,
	xlab = 'True Simulated Values', ylab = 'Predicted Values',
	cex.lab = 2, cex.axis = 1.5)
# make a map (bubble plot) of the predicted values
par(mar = c(5,5,1,1))
plot(p1$sx, p1$sy, cex = predfit$fit,
	xlab = 'x-coordinate', ylab = 'y-coordinate')
points(d1$sx, d1$sy, cex = d1$y, pch = 19)

# ------------------------------------------------------------------------------
# Model Selection 
# ------------------------------------------------------------------------------

# create a fixed spatial pattern
# first, create 10,000 spatial locations on a grid
sx = as.vector(outer(1:100, rep(1, times = 100)))/100
sy = as.vector(outer(rep(1, times = 100), 1:100))/100
par(mar = c(5,5,1,1))
plot(sx, sy, pch = 19, cex = .3, cex.lab = 2, cex.axis = 1.5)
# do a rotation on the spatial coordinates
rot = .1
rotc = cbind(sx, sy) %*% matrix(c(cos(rot*pi), -sin(rot*pi), 
	-sin(rot*pi), cos(rot*pi)), nrow = 2, byrow = TRUE)
plot(rotc)
# create a 2-D sine wave surface
z = exp(.4*sin(0.4*2*pi*(20*rotc[,1] + 0.9*pi)) + 
	0.4*sin(0.1*2*pi*(20*rotc[,2] + 0.9*pi)))
z = z/max(z)

nclass = 12
cip = classIntervals(z, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)
# add a bunch of these sine waves together at various rotations,
# frequencies, and amplitudes
z = rep(0, times = 100^2)
for(i in 1:20) {
	rot = i/21
	rotc = cbind(sx, sy) %*% matrix(c(cos(rot*pi), -sin(rot*pi), 
		-sin(rot*pi), cos(rot*pi)), nrow = 2, byrow = TRUE)
	# create a 2-D sine wave surface
	z = z + exp((i/21)*sin(((21 - i)/21)*2*pi*(20*rotc[,1] + (i/21)*pi)) + 
		((21 - i)/21)*sin((i/21)*2*pi*(20*rotc[,2] + ((21 - i)/21)*pi)))
}
z = z/(max(z))
nclass = 12
cip = classIntervals(z, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# create some smooth bell-shaped surfaces
z1 = exp(-((sx - 0.3)^2 + (sy - 0.3)^2)/.2)
nclass = 12
cip = classIntervals(z1, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

z2 = .6*exp(-((sx - 0.1)^2 + (sy - 0.9)^2)/.05)
nclass = 12
cip = classIntervals(z2, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# add the two bell-shaped surfaces together
z12 = z1 + z2
cip = classIntervals(z12, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# add sine-wave and bell-shaped surfaces together
z = z + z1 + z2
# standardize
z = (z - mean(z))/sd(z)
cip = classIntervals(z, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# create 2 sets of x-variables, assume fixed by using set.seed()
set.seed(3009)
x1 = rnorm(100^2)
x3 = rnorm(100^2)

# create a nonlinear response to x1 using sines and exponentials
y1 = .1*sin(x1/.02) + .2*sin(x1/.05) + .3*sin(x1/.1) + 
	exp(sin(x1/3))
y1 = (y1 - mean(y1))/sd(y1)
par(mar = c(5,5,1,1))
plot(x1, y1, pch = 19, cex.lab = 2, cex.axis = 1.5)

# create x2 as a spatial surface
x2 = .1*sin(sx/.02) + .2*sin(sx/.05) + .3*sin(sx/.1) + 
	exp(sin(sx/3)) + .1*sin(sy/.01) + .2*sin(sy/.01) + .3*sin(sy/.05) + 
	exp(sin(sy/2)) 
x2 = (x2 - mean(x2))/sd(x2)
cip = classIntervals(x2, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# create our response from y1, x2, and x3
y = y1 + 0*x2 + 0.1*x3 + z

# plot the response variable in space
cip = classIntervals(y, n = nclass, style = 'fisher')
palp = viridis(nclass)
cont_colors = findColours(cip, palp)
plot(sx, sy, pch = 19, col = cont_colors)

# create data and prediction data sets
set.seed(4010)
samp = sample(1:10000, size = 1100)
d1 = data.frame(sx = sx[samp[1:100]], sy = sy[samp[1:100]],
	x1 = x1[samp[1:100]], x2 = x2[samp[1:100]], 
	y = y[samp[1:100]])
p1 = data.frame(sx = sx[samp[101:1100]], sy = sy[samp[101:1100]],
	x1 = x1[samp[101:1100]], x2 = x2[samp[101:1100]], 
	y = y[samp[101:1100]])

# suppose that you can only measure x1 and x2
# scatterplot of y versus x1
par(mar =c(5,5,1,1))
plot(d1$x1, d1$y, pch = 19, cex.lab = 2, cex.axis = 1.5, cex = 2)
# scatterplot of y versus x2
par(mar =c(5,5,1,1))
plot(d1$x2, d1$y, pch = 19, cex.lab = 2, cex.axis = 1.5, cex = 2)

# lm and splm give the same results, except splm assumes z-distribution
lmout = lm(y ~ x1 + x2, data = d1)
summary(lmout)
splmout = splm(y ~ x1 + x2, data = d1, spcov_type = 'none')
summary(splmout)

# fit several spatial models with and without x2
splmexp = splm(y ~ x1 + x2, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'exponential', estmethod = 'ml')
summary(splmexp)
splmsph = splm(y ~ x1 + x2, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'spherical', estmethod = 'ml')
summary(splmsph)
splmgau = splm(y ~ x1 + x2, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'gaussian', estmethod = 'ml')
summary(splmgau)
splmexp1 = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'exponential', estmethod = 'ml')
summary(splmexp1)
splmsph1 = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'spherical', estmethod = 'ml')
summary(splmsph1)
splmgau1 = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'gaussian', estmethod = 'ml')
summary(splmgau1)

# compare models using AIC
AIC(splmout)
AIC(splmexp)
AIC(splmsph)
AIC(splmgau)
AIC(splmexp1)
AIC(splmsph1)
AIC(splmgau1)

# use glances to make it easy, quick, and offer other comparisons
glances(splmout, splmexp, splmsph, splmgau, splmexp1, splmsph1, splmgau1)

# compare all of the models using LOOCV
loocv(splmout)
loocv(splmexp)
loocv(splmsph)
loocv(splmgau)
loocv(splmexp1)
loocv(splmsph1)
loocv(splmgau1)

# make final fits of the models using REML rather than ML
splmexp1reml = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'exponential', estmethod = 'reml')
summary(splmexp1reml)
splmsph1reml = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'spherical', estmethod = 'reml')
summary(splmsph1reml)
splmgau1reml = splm(y ~ x1, data = d1, xcoord = sx, ycoord = sy,
	spcov_type = 'gaussian', estmethod = 'reml')
summary(splmgau1reml)

# what are the models trying to estimate with linear model?
popfit = lm(y1 ~ x1)$coefficient
popfit[2]
plot(x1, y1, pch = 19, cex.lab = 2, cex.axis = 1.5)
lines(c(-4, 4), c(popfit[1] - 4*popfit[2], popfit[1] + 4*popfit[2]),
	lwd = 3)

# compare models after using REML
glances(splmexp1reml, splmsph1reml, splmgau1reml)

# make predictions for the unobserved data
indpreds = predict(splmout, p1)
exppreds = predict(splmexp1reml, p1)
sphpreds = predict(splmsph1reml, p1)
gaupreds = predict(splmgau1reml, p1)

# compare the predictions to the true values
par(mar = c(5,5,1,1))
plot(p1$y, indpreds, pch = 19, cex.lab = 2, cex.axis = 1.5)
plot(p1$y, exppreds, pch = 19, cex.lab = 2, cex.axis = 1.5)

# compare RMSPE for predictions versus true values
sqrt(mean((indpreds - p1$y)^2))
sqrt(mean((exppreds - p1$y)^2))
sqrt(mean((sphpreds - p1$y)^2))
sqrt(mean((gaupreds - p1$y)^2))
