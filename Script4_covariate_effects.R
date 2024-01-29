# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Bayesian band-recovery workshop
# North American Duck Symposium
#
# covariate example: S(t) f(Duck stamps)
#                    S(t) f(Duck stamps + t)
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
# Duck Stamps (DS)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- NULL
D[1] <- 0
for (t in 2:(nT)){
  # this is an 'autoregressive' parameterization that will lead to several 
  # years of higher hunter numbers followed by several years of low hunter numbers
  D[t] <- rnorm(1, D[t-1] * 0.5, 0.25)
}
plot(D, type = 'b', ylab = 'Duck Stamps')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival probability for each year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
muS <- qlogis(0.6) # 0.5 is the mean S on the prob. scale
sdS <- 0.25        # 0.25 is the st dev of S on the logit scale
epsilon.S <- rnorm(nT-1, 0, sdS)
S <- plogis(muS + epsilon.S)

plot(S ~ seq(1,nT-1), type = 'b', xlab = 'T')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# band-recovery probability (f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
muf <- qlogis(0.05)  # 0.1 is the mean f on the prob. scale
sdf <- 0.1           # 0.15 is the st dev of f on the logit scale
beta <- 1
epsilon.f <- rnorm(nT, 0, sdf)
f <- plogis(muf + beta * D + epsilon.f)

plot(f ~ seq(1,nT), type = 'b', xlab = 'T', ylim = c(0,0.15))
plot(f ~ D)


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
# S(t) f(Duck stamps)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m4a.jags")
cat("
      model {

      for (t in 1:(nT-1)){
        S[t] ~ dbeta(1,1)
      }
      
      alpha ~ dlogis(0,1)
      beta ~ dnorm(0,0.1)
      for (t in 1:nT){
        logit(f[t]) = alpha + beta * D[t]
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
jags.data <- list(marr = marr, rel = rowSums(marr), nT = nT, D = D)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f','alpha','beta')

nc <- 4
nt <- 10
ni <- 15000
nb <- 10000


Sys.time()
m4a <- jags(jags.data, inits, parameters, "m4a.jags", 
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
plot(m4a$q50$f ~ D, 
     ylim = c(0, 0.15), 
     ylab = 'Band recovery probability (f)',
     xlab = 'Duck Stamps')
arrows(D, m4a$q2.5$f, D, m4a$q97.5$f,
       lty = 2, length = 0)
points(m4a$q50$f, pch = 21, bg = 'white', cex = 1.25)
points(f ~ D, pch = 19, col = 'red')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot intercept, covariate effect, and proportion of posterior distribution
# less than 0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vioplot(plogis(m4a$sims.list$alpha))
vioplot(m4a$sims.list$beta)

# proportion of posterior distribution less than 0
m4a$f$beta
length(which(m4a$sims.list$beta < 0))/length(m4a$sims.list$beta)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting S estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m4a$q50$S, 
     ylim = c(0.3,0.9), 
     ylab = 'Survival probability (f)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT), m4a$q2.5$S, 1:(nT), m4a$q97.5$S,
       lty = 2, length = 0)
points(m4a$q50$S, pch = 21, bg = 'grey50', cex = 1.25)
points(S, pch = 19, type = 'b', col = 'red')












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# S(t) f(Duck stamps + t)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m4b.jags")
cat("
      model {

      for (t in 1:(nT-1)){
        S[t] ~ dbeta(1,1)
      }
      
      alpha ~ dlogis(0,1)
      beta ~ dnorm(0,0.1)
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)
      
      for (t in 1:nT){
        epsilon[t] ~ dnorm(0, tau)
        logit(f[t]) = alpha + beta * D[t] + epsilon[t]
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
jags.data <- list(marr = marr, rel = rowSums(marr), nT = nT, D = D)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f','alpha','beta','epsilon')

nc <- 4
nt <- 10
ni <- 15000
nb <- 10000


Sys.time()
m4b <- jags(jags.data, inits, parameters, "m4b.jags", 
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
cov <- seq(-1, 1, length.out = 100)
pred <- plogis(m4b$q50$alpha + m4b$q50$beta * cov)

plot(m4b$q50$f ~ D, 
     ylim = c(0, 0.15), 
     ylab = 'Band-recovery probability (f)',
     xlab = 'Duck stamps (D)')
arrows(D, m4b$q2.5$f, D, m4b$q97.5$f,
       lty = 2, length = 0)
points(m4b$q50$f ~ D, pch = 21, bg = 'white', cex = 1.25)
points(f ~ D, pch = 19, col = 'red')
points(pred ~ cov, type = 'l', col = 'purple', lty = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot intercept, covariate effect, and proportion of posterior distribution
# less than 0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vioplot(plogis(m4b$sims.list$alpha))
vioplot(m4b$sims.list$beta)

# proportion of posterior distribution less than 0
m4b$f$beta
length(which(m4b$sims.list$beta < 0))/length(m4b$sims.list$beta)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting band-recovery estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m4b$q50$S, 
     ylim = c(0,1), 
     ylab = 'Survival probability (S)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m4b$q2.5$S, 1:(nT-1), m4b$q97.5$S,
       lty = 2, length = 0)
points(m4b$q50$S, pch = 21, bg = 'grey50', cex = 1.25)
points(S, pch = 19, type = 'b', col = 'red')




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 1) Our estimates of the effects of duck stamps on band-recovery probability vary between the two model
#    types (see the plot). How and Why? 
vioplot(m4a$sims.list$beta, m4b$sims.list$beta)
#
# 2) Why would we wish to estimate covariate effects as well as additional
#    random temporal variation simultaneously? Try to run a model without the covariate (random-effects only)
#    What happens to the estimate of the variance for f
# 
# 3) What is the effect of changing sample sizes on the precision of our estimates?
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




