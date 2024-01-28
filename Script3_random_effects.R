# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Bayesian band-recovery workshop
# North American Duck Symposium
#
# random effects example: S(t) f(t)
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

pr[1,2]
S[1] * f[2]

pr[3,5]
S[3] * S[4] * f[5]

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
sink("m3.jags")
cat("
      model {

      sigma.f ~ dgamma(1,1)
      sigma.S ~ dgamma(1,1)

      tau.f = 1/(sigma.f * sigma.f)
      tau.S = 1/(sigma.S * sigma.S)

      mean.f ~ dlogis(0,1)
      mean.S ~ dlogis(0,1)

      for (t in 1:nT){
        epsilon.f[t] ~ dnorm(0, tau.f)
        logit(f[t]) = mean.f + epsilon.f[t]
      }
      
      for (t in 1:(nT-1)){
        epsilon.S[t] ~ dnorm(0, tau.S)
        logit(S[t]) = mean.S + epsilon.S[t]
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
parameters <- c('S','f',
                'mean.f','mean.S',
                'sigma.f','sigma.S',
                'epsilon.f','epsilon.S')

nc <- 4
nt <- 10
ni <- 15000
nb <- 10000


Sys.time()
m3 <- jags(jags.data, inits, parameters, "m3.jags", 
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
plot(m3$q50$S, 
     ylim = c(0,1), 
     ylab = 'Survival (S)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m3$q2.5$S, 1:(nT-1), m3$q97.5$S,
       lty = 2, length = 0)
points(m3$q50$S, pch = 21, bg = 'white', cex = 1.25)
points(S, pch = 19, type = 'b', col = 'red')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting band-recovery estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m3$q50$f, 
     ylim = c(0,0.25), 
     ylab = 'Band-recovery probability (f)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m3$q2.5$f, 1:(nT-1), m3$q97.5$f,
       lty = 2, length = 0)
points(m3$q50$f, pch = 21, bg = 'grey50', cex = 1.25)
points(f, pch = 19, type = 'b', col = 'red')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 1) As we discussed, using random effects can lead to 'shrinkage'
#    of annual estimates towards the mean in years where sample sizes
#    are small, or when sample sizes are small over the entire study.
#
#    a) can you think of a situation in which this might be desireable?
#    b) would this ever be undesireable? why?
#
#  2) In your own words, why are they called 'random' effects? What
#     biological hypothesis does this represent? When might this be
#     a poor hypothesis?
#
#  3) Many researchers recommend a minimum sample size of 7 (or 10)
#     or more 'groups' (or years) when using random effects. Why?
#     Perhaps consider a special case of 2 years of data as an extreme
#     violation of this rule? How might you estimate the properties of a
#     distribution with two samples?
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

