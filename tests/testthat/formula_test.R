rm(list = ls())
set.seed(12345)
n <- 100
p <- 4
Y <- rnorm(n)
d <- rnorm(n)
X <- matrix(rnorm(n * p), ncol = p)
Z <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste("X", 1:p, sep = "")
colnames(Z) <- paste("Z", 1:p, sep = "")
mat <- cbind(Y, d, X, Z)
colnames(mat)[1:2] <- c("y", "d")
dat <- as.data.frame(mat)

form <- as.Formula(Y ~ d + X1 + X2 + X3 + X4 | Z1 + Z2 + Z3 + Z4 + X1 + X2 + X3 + X4)
dataf <- f.formula(Y ~ d + X1 + X2 + X3 + X4 | X1 + X2 + X3 + X4, data = dat)
dataf <- f.formula(Y ~ d + X1 + X2 + X3 + X4 | Z1 + Z2 + Z3 + Z4 + X1 + X2 + X3 + X4, data = dat)
dataf$Y
dataf$D
head(dataf$X)
head(dataf$Z)

## test function include
i <- include(I = ~ X1 + X2, dat)
test1 <- rlassoEffects(Y ~ X1 + X2 + Z1 + Z2, I = ~X1, data = dat)
# test2 <- rlassoEffects(Y ~ X1 + X2+ Z1 + Z2, I= ~ X1)
test3 <- rlassoEffects(y ~ X1 + X2 + Z1 + Z2, I = ~X1, include = ~X1, data = dat)
test4 <- rlassoEffects(y ~ X1 + X2 + Z1 + Z2, I = ~X1, include = ~X3, data = dat)
test5 <- rlassoEffects(y ~ (X1 + X2)^2 + Z1 + Z2, I = ~X1, include = ~X3, data = dat)
test6 <- rlassoEffects(y ~ (X1 + X2)^2 + Z1 + Z2, I = ~ X1:X2, include = ~X3, data = dat)
test7 <- rlassoEffects(y ~ (X1 + X2)^2 + Z1 + Z2, I = ~ X1 + X2, include = ~X3, data = dat)
test8 <- rlassoEffects(y ~ (X1 + X2)^2 + Z1 + Z2, I = ~ X1 + X2, include = ~NULL, data = dat)




############################ Test IV functions
test1a <- rlassoIVselectX(form, data = dat)
test1b <- rlassoIVselectX(x = X, z = Z, y = Y, d = d)
test2a <- rlassoIVselectZ(form, data = dat)
test2b <- rlassoIVselectZ(x = X, z = Z, y = Y, d = d)
# testold <- rlassoIVselectZ(x=X, d=d, y=Y, z=Z)
test3a <- rlassoIV(form, data = dat)
test3b <- rlassoIV(x = X, z = Z, y = Y, d = d)
test4a <- tsls(form, data = dat)
test4b <- tsls(x = X, z = Z, y = Y, d = d)

all.equal(test1a, test1b)
all.equal(test2a, test2b)
all.equal(test3a, test3b)
all.equal(test4a, test4b)

############################################################################
## now d multivariate!!!!!!!!!!!!!
############################################################################

set.seed(12345)
n <- 100
p <- 5
Y <- rnorm(n)
d1 <- rnorm(n)
d2 <- rnorm(n)
d <- cbind(d1, d2)
colnames(d) <- c("d1", "d2")
X <- matrix(rnorm(n * p), ncol = p)
Z <- matrix(rnorm(n * p), ncol = p)
colnames(X) <- paste("X", 1:p, sep = "")
colnames(Z) <- paste("Z", 1:p, sep = "")
mat <- cbind(Y, d1, d2, X, Z)
colnames(mat)[1:3] <- c("y", "d1", "d2")
dat <- as.data.frame(mat)

form <- as.Formula(Y ~ d1 + d2 + X1 + X2 + X3 + X4 | Z1 + Z2 + Z3 + Z4 + X1 + X2 + X3 + X4)


#########################################################################################
dataf <- f.formula(Y ~ d1 + d2 + X1 + X2 + X3 + X4 | X1 + X2 + X3 + X4, data = dat)
dataf <- f.formula(Y ~ d1 + d2 + X1 + X2 + X3 + X4 | Z1 + Z2 + Z3 + Z4 + X1 + X2 + X3 + X4, data = dat)
dataf$Y
dataf$D
head(dataf$X)
head(dataf$Z)

############################ Test IV functions
test1a <- rlassoIVselectX(form, data = dat)
test1b <- rlassoIVselectX(x = X, z = Z, y = Y, d = d) # dim(d)[2]==1 !!!
test2a <- rlassoIVselectZ(form, data = dat)
test2b <- rlassoIVselectZ(x = X, z = Z, y = Y, d = d) # dim(d)[2]>1 works in principle
test3a <- rlassoIV(form, data = dat)
test3b <- rlassoIV(x = X, z = Z, y = Y, d = d) # dim(d)[2]==1 !!!
test4a <- tsls(form, data = dat)
test4b <- tsls(x = X, z = Z, y = Y, d = d) # dim(d)[2]>1 works in principle

all.equal(test1a, test1b)
all.equal(test2a, test2b)
all.equal(test3a, test3b)
all.equal(test4a, test4b)


#############################################################################################################
############################### Test treatment effects
#############################################################################################################


set.seed(12345)
n <- 100
p <- 4
Y <- rnorm(n)
d <- sample(c(0, 1), size = n, replace = TRUE)
X <- matrix(rnorm(n * p), ncol = p)
Z <- sample(c(0, 1), size = n, replace = TRUE)
colnames(X) <- paste("X", 1:p, sep = "")
mat <- cbind(Y, d, Z, X)
colnames(mat)[1:3] <- c("y", "d", "Z")
dat <- as.data.frame(mat)

form <- as.Formula(Y ~ d + X1 + X2 + X3 + X4 | Z + X1 + X2 + X3 + X4)
form2 <- as.Formula(Y ~ d + X1 + X2 + X3 + X4 | X1 + X2 + X3 + X4)
test1a <- rlassoATE(form2, data = dat)
test1b <- rlassoATE(x = X, y = Y, d = d)
test2a <- rlassoATET(form2, data = dat)
test2b <- rlassoATET(x = X, y = Y, d = d)
test3a <- rlassoLATE(form, data = dat)
test3b <- rlassoLATE(x = X, z = Z, y = Y, d = d)
test4a <- rlassoLATET(form, data = dat)
test4b <- rlassoLATET(x = X, z = Z, y = Y, d = d)

all.equal(test1a, test1b)
all.equal(test2a, test2b)
all.equal(test3a, test3b)
all.equal(test4a, test4b)












############################### Test treatment effects

form2 <- as.Formula(Y ~ d1 + d2 + X1 + X2 + X3 + X4 | X1 + X2 + X3 + X4)
test1a <- rlassoATE(form2, data = dat)
test1b <- rlassoATE(x = X, y = Y, d = d)
test2a <- rlassoATET(form2, data = dat)
test2b <- rlassoATET(x = X, y = Y, d = d)
test3a <- rlassoLATE(form, data = dat)
test3b <- rlassoLATE(x = X, z = Z, y = Y, d = d)
test4a <- rlassoLATET(form, data = dat)
test4b <- rlassoLATET(x = X, z = Z, y = Y, d = d)

all.equal(test1a, test1b)
all.equal(test2a, test2b)
all.equal(test3a, test3b)
all.equal(test4a, test4b)



##################################################################################################################################
