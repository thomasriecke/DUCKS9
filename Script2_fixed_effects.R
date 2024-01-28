# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Bayesian band-recovery workshop
# North American Duck Symposium
#
# fixed effects example: S(t) f(t)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
library(vioplot)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'serif')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of years (nT)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nT <- 15

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of releases each year (nR)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nR <- 1000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival probability for each year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
muS <- qlogis(0.5) # 0.5 is the mean S on the prob. scale
sdS <- 0.25        # 0.25 is the st dev of S on the logit scale

lS <- rnorm(nT-1, muS, sdS)
S <- plogis(lS)

plot(S ~ seq(1,nT-1), type = 'b', xlab = 'T', ylim = c(0.3,0.7))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# band-recovery probability (f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
muf <- qlogis(0.1) # 0.1 is the mean f on the prob. scale
sdf <- 0.15        # 0.15 is the st dev of f on the logit scale

lf <- rnorm(nT, muf, sdf)
f <- plogis(lf)

plot(f ~ seq(1,nT), type = 'b', xlab = 'T', ylim = c(0,0.25))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualizing the logit-link
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(plogis(seq(-3,3,length.out = 100)) ~ seq(-3,3,length.out = 100),
     type = 'l', ylab = 'S', xlab = 'logit(S)', cex.lab = 1.5)
points(f ~ lf)
points(S ~ lS, pch = 19)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define the cell probabilities (pr) for the m-array
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pr <- matrix(NA, nT, nT+1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fill the main diagonal (direct recovery probability)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (t in 1:nT){
  pr[t,t] <- f[t]
}
f
View(pr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fill below the main diagonal (all 0's!)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (t in 2:nT){
  for (j in 1:(t-1)){
    pr[t,j] <- 0
  }
}
View(pr)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fill above the main diagonal
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (t in 1:(nT-1)){
  for (j in (t+1):nT){
    pr[t,j] <- prod(S[t:(j-1)]) * f[j]
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# last column
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (t in 1:nT){
  pr[t,(nT+1)] <- 1 - sum(pr[t,1:nT])
}
rowSums(pr)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate band-recovery data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
marr <- matrix(NA, nT, nT+1)
for (t in 1:nT){
  marr[t,] <- rmultinom(1, nR, pr[t,])
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# existing model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m2.jags")
cat("
      model {

      for (t in 1:nT){
        f[t] ~ dbeta(1,1)
      }
      
      for (t in 1:(nT-1)){
        S[t] ~ dbeta(1,1)
      }
      


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # m-array likelihood
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (t in 1:nT){
      marr[t,] ~ dmulti(pr[t,], rel[t])
      pr[t,t] = f[t]
      pr[t,(nT+1)] <- 1 - sum(pr[t,1:nT])
    }

    for (t in 1:(nT-1)){
      for (j in (t+1):nT){
        pr[t,j] = prod(S[t:(j-1)]) * f[j]
      }
    }
    
    for (t in 2:nT){
      for (j in 1:(t-1)){
        pr[t,j] = 0
      }
    }

      }
      ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(marr = marr, rel = rowSums(marr), nT = nT)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f')

nc <- 4
nt <- 10
ni <- 15000
nb <- 10000


Sys.time()
m2 <- jags(jags.data, inits, parameters, "m2.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot medians and confidence intervals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting survival estimates (white) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m2$q50$S, 
     ylim = c(0,1), 
     ylab = 'Survival (S)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m2$q2.5$S, 1:(nT-1), m2$q97.5$S,
       lty = 2, length = 0)
points(m2$q50$S, pch = 21, bg = 'white', cex = 1.25)
points(S, pch = 19, type = 'b', col = 'red')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting band-recovery estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m2$q50$f, 
     ylim = c(0,0.25), 
     ylab = 'Band-recovery probability (f)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m2$q2.5$f, 1:(nT-1), m2$q97.5$f,
       lty = 2, length = 0)
points(m2$q50$f, pch = 21, bg = 'grey50', cex = 1.25)
points(f, pch = 19, type = 'b', col = 'red')








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions
# 
# 1) Why don't the estimates perfectly match the true underlying parameter
#    estimates? How might this change if we changed the number
#    of released individuals (up or down)? Why?
#
# 2) Why is there one more band-recovery probability estimate (nT)
#    than survival estimate (nT-1)?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


