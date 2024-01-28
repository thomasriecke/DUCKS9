# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bayesian band-recovery workshop
# Duck Symposium
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
library(vioplot)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up the simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nT <- 10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of releases each year
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nR <- 1000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# hunting mortality probability (kappa)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kappa <- plogis(rnorm(nT, -1.5, 0.5))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# crippling loss probability (c)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c <- 0.2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# harvest mortality probability (h)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h <- kappa * (1 - c)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# band-reporting probability (rho)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rho <- 0.5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# band-recovery probability (f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f <- h * rho

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define the cell probabilities (pr) for the m-array
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r <- rbinom(nT, nR, f)


sink("hunting_harvest.jags")
cat("
    model {

    rho ~ dbeta(100, 100)
    c ~ dbeta(20,80)
    
    for (t in 1:nT){
      kappa[t] ~ dbeta(1,1)
      h[t] = kappa[t] * (1 - c)
      f[t] = kappa[t] * rho * (1 - c)
    }
 
 
    for (t in 1:nT){
      r[t] ~ dbin(f[t], nR)
    } 

      }
      ",fill = TRUE)
sink()

#############################################################################################
# Bundle data
#############################################################################################

jags.data <- list(r = r, nR = nR, nT = nT)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('f','rho','h', 'kappa','c')

nc <- 4
nt <- 10
ni <- 15000
nb <- 10000

library(jagsUI)
Sys.time()
m1 <- jags(jags.data, inits, parameters, "hunting_harvest.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()

print(m1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting band-recovery estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m1$q50$f, 
     ylim = c(0,0.25), 
     ylab = 'Band-recovery probability (f)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m1$q2.5$f, 1:(nT-1), m1$q97.5$f,
       lty = 2, length = 0)
points(m1$q50$f, pch = 21, bg = 'grey50', cex = 1.25)
points(f, pch = 19, type = 'b', col = 'red')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting harvest estimates (grey) against 'truth' (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(m1$q50$h, 
     ylim = c(0,0.4), 
     ylab = 'Harvest probability (f)',
     xlab = 'Occasion (T)', 
     type = 'b')
arrows(1:(nT-1), m1$q2.5$h, 1:(nT-1), m1$q97.5$h,
       lty = 2, length = 0)
points(m1$q50$h, pch = 21, bg = 'grey50', cex = 1.25)
points(h, pch = 19, type = 'b', col = 'red')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# violin plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vioplot(m1$sims.list$rho, m1$sims.list$c,
        drawRect = F, wex = 0.5,
        ylab = 'Reporting and crippling rate',
        names = c(expression(rho),'c'), 
        ylim = c(0,0.65), cex.lab = 2)
arrows(1:2, c(m1$q2.5$rho,m1$q2.5$c), 
       1:2, c(m1$q97.5$rho,m1$q97.5$c), length = 0, col = 'white', lwd = 2)
points(c(m1$q50$rho, m1$q50$c), cex = 2, pch = 21, bg = 'white')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions:
# 1) make a plot of hunting mortality probability (kappa)
#    similar to the plots for
#    band-recovery probability (f) and harvest probability (h)
#
# 2) Change the prior for band-reporting rate to dbeta(40,60)
#    what will this do to our estimate of band-reporting probability?
#    how will this affect our estimates of kappa?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
