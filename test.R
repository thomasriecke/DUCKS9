library(jagsUI)
library(vioplot)

y <- 10
N <- 100

sink("test.jags")
cat("
      model {

       theta ~ dbeta(1,1)
       y ~ dbin(theta, N)
   
      }
      ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bundle data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(y = y, N = N)

# Initial values
inits <- function(){list(theta = 0.1)}  

# Parameters monitored
parameters <- c('theta')

nc <- 4
nt <- 1
ni <- 5000
nb <- 1000


Sys.time()
library(jagsUI)
m <- jags(jags.data, inits, parameters, "test.jags", 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
          parallel = T)
Sys.time()
print(m, digits = 3)

vioplot(m$sims.list$theta, drawRect = F, col = 'royalblue1', names = '')
arrows(1, m$q2.5$theta, 1, m$q97.5$theta, length = 0, col = 'white')
arrows(1, m$q25$theta, 1, m$q75$theta, length = 0, col = 'white', lwd = 5)
points(m$q50$theta ~ 1, cex = 3, pch = 21, bg = 'white')
mtext(expression(theta), side = 2, line = 3, cex = 1.5)
