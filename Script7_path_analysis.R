# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this script reads AHY female midcontinent mallard banding data from 1974-1994,
# and then generates estimates of abundance using harvest estimates
#
# D: Duck stamp sales
# N: BPOP abundance estimates (in thousands)
# P: Ponds in the US and Canada (in thousands)
# 
# 
# ~ Thomas Riecke 
# thomasvanceriecke@gmail.com
# thomas.riecke@umontana.edu
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x <- url("https://github.com/thomasriecke/DUCKS9/blob/main/mallard_data.RData?raw=true")
load(x)


# z-standardizing the BPOP data and Duck stamp sales to help with convergence.
zN <- (N - mean(N))/sd(N)
zD <- (D - mean(D))/sd(D)
plot(zD ~ zN)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-seasonal model with Lincoln estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(jagsUI)
# Specify model in BUGS language
sink("path_analysis_hunting.jags")
cat("
    model {




    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    #
    # band-reporting rate and crippling loss priors
    #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    c ~ dbeta(20,80)
    for (t in 1:n.years){
      rho[t] ~ dbeta(rho.alpha[t], rho.beta[t])
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # random-effects parameterization
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alpha.phi ~ dlogis(0,1)
    alpha.omega ~ dlogis(0,1)
    alpha.kappa ~ dlogis(0,1)
    alpha.D ~ dnorm(0,1)
    
    sigma.D ~ dgamma(0.1, 0.1)
    sigma.phi ~ dgamma(1,1)
    sigma.omega ~ dgamma(1,1)
    sigma.kappa ~ dgamma(1,1)
    
    tau.D = 1/(sigma.D * sigma.D)
    tau.phi = 1/(sigma.phi * sigma.phi)
    tau.omega = 1/(sigma.omega * sigma.omega)
    tau.kappa = 1/(sigma.kappa * sigma.kappa)

    beta.kappa[1] ~ dnorm(0, 0.1)
    beta.kappa[2] ~ dnorm(0, 0.1)
    beta.D ~ dnorm(0, 0.1)


    for(t in 1:n.years){
    
      epsilon.phi[t] ~ dnorm(0, tau.phi)
      logit(phi[t]) = alpha.phi + epsilon.phi[t] 
      
      epsilon.kappa[t] ~ dnorm(0, tau.kappa)
      logit(kappa[t]) = alpha.kappa + beta.kappa[1] * zN[t] + 
                        beta.kappa[2] * zD[t] + epsilon.kappa[t]
      
      zD[t] ~ dnorm(alpha.D + beta.D * zN[t], tau.D)
      
      f[t] = kappa[t] * (1 - c) * rho[t]
      h[t] = kappa[t] * (1 - c)
      
      
    }

    for (t in 1:(n.years-1)){

      epsilon.omega[t] ~ dnorm(0, tau.omega)
      logit(omega[t]) = alpha.omega + epsilon.omega[t]
      S[t] = phi[t+1] * (omega[t] * (1 - kappa[t]))
      
    }




    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # adult female post-season m-array [,,1]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (t in 1:n.years){
        marr[t, 1:(n.years+1), 1] ~ dmulti(pr[t,,1], rel[t,1])
        pr[t,n.years+1,1] <- 1-sum(pr[t, 1:n.years, 1])
        pr[t,t,1] <- phi[t] * f[t] 
    }

    for (t in 1:(n.years-1)){
      for (j in (t+1):n.years){
        pr[t,j,1] <- phi[t] * prod(S[t:(j-1)]) * f[j]
      }
    }

    for (t in 2:n.years){
      for (j in 1:(t-1)){
        pr[t,j,1] <- 0 
      }
    }
    

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # adult female pre-season m-array [,,2]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (t in 1:n.years){
        marr[t, 1:(n.years+1), 2] ~ dmulti(pr[t,,2], rel[t,2])
        pr[t,n.years+1,2] <- 1-sum(pr[t,1:n.years,2])
        pr[t,t,2] <- f[t] 
    }

    for (t in 1:(n.years-1)){
      for (j in (t+1):n.years){
        pr[t,j,2] <- prod(S[t:(j-1)]) * f[j]
      }
    }

    for (t in 2:n.years){
      for (j in 1:(t-1)){
        pr[t,j,2] <- 0 
      }
    }


    

    }
    ",fill = TRUE)
sink()

#############################################################################################
# Bundle data
#############################################################################################
jags.data <- list(marr = marr, rel = rel, n.years = n.years, 
                  zN = zN, zD = zD,
                  rho.alpha = rho.alpha, rho.beta = rho.beta)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f','h','phi','omega','rho','kappa','c','alpha.D','alpha.kappa',
                'sigma.h','sigma.omega','sigma.phi','beta.D','beta.kappa',
                'epsilon.h','epsilon.omega','epsilon.phi','sigma.D')

nc <- 4
nt <- 25
ni <- 10000
nb <- 2500


# Call JAGS from R (BRT )
Sys.time()
m7 <- jags(jags.data, inits, parameters, "path_analysis_hunting.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot estimates of hunting mortality against abundance (N) and
# duck stamp sales (D)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
greens <- colorRampPalette(c('lightgreen','darkgreen'))
plot(m7$q50$kappa ~ N, ylab = expression(Hunting~mortality~probability~(kappa)),
     xlab = 'WBPHS N (millions)', pch = 21, bg = greens(n.years), cex = 2)
plot(m7$q50$kappa ~ D, ylab = expression(Hunting~mortality~probability~(kappa)),
     xlab = 'Duck stamp sales', pch = 21, bg = greens(n.years), cex = 2)

vioplot(m7$sims.list$beta.kappa[,1], m7$sims.list$beta.kappa[,2], m7$sims.list$beta.D,
        drawRect = FALSE, col = 'forestgreen',
        names = c(expression(beta[kappa*'1']),expression(beta[kappa*'2']),expression(beta[D])))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a figure demonstrating relationship between kappa and 
# duck stamp sales
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a range of covariate values (-2, 2)
res <- 100
x <- seq(-2,2,length.out = res)

# generate predicted values given uncertainty in posterior distributions
pred.kappa <- matrix(NA, length(m7$sims.list$alpha.kappa), res)
for (j in 1:res){
  pred.kappa[,j] <- plogis(m7$sims.list$alpha.kappa + m7$sims.list$beta.kappa[,2] * x[j])
}

# generate 90% quantiles of predictions
q.pred <- matrix(NA, res, 3)
for (j in 1:res){
  q.pred[j,] <- quantile(pred.kappa[,j], c(0.05,0.5,0.95))
}

# melt the predictions
library(reshape2)
p.kap <- melt(pred.kappa)
names(p.kap) <- c('index','x','est')
x <- seq(-2,2,length.out = 5)
x.labs <- x * sd(D) + mean(D)

# plot
greens <- colorRampPalette(c('white','forestgreen','darkgreen'))
smoothScatter(p.kap$est ~ p.kap$x, nrpoints = 0, 
              ylab = expression(Hunting~mortality~probability~(kappa)),
              xlab = 'Duck stamp sales', colramp = greens,
              xaxt = 'n', las = 1, cex.lab = 1.5)

axis(side = 1, at = seq(0,100,length.out = 5), 
     labels = seq(min(D), max(D), length.out = 5))
lines(q.pred[,2] ~ seq(1,100), col = 'white', lwd = 3)
lines(q.pred[,1] ~ seq(1,100), col = 'white', lwd = 2, lty = 3)
lines(q.pred[,3] ~ seq(1,100), col = 'white', lwd = 2, lty = 3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a similar figure using estimates and raw data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
greens <- colorRampPalette(c('lightgreen','forestgreen','darkgreen'))
plot(m7$q50$kappa ~ D, ylab = expression(Hunting~mortality~probability~(kappa)),
     xlab = 'Duck stamp sales', pch = 21, cex = 1, cex.lab = 1.5,
     ylim = c(0,0.15))
arrows(D, m7$q2.5$kappa, D, m7$q97.5$kappa, length = 0, lty = 2, col = greens(n.years))
points(m7$q50$kappa ~ D, cex = 2, bg = greens(n.years), pch = 21)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions
# 1) If we simply plot estimates of kappa (or f or h) against N, we see
#    a strong positive relationship. What does our hypothesized path
#    path analysis indicate about this relationship? Is this a direct effect
#    of changes in N, or is this mediated via another causal mechanism?
#
# 2) How might our management actions change in light of this information.
#    In other words. Imagine a previous analysis that had found a strong linkage
#    between kappa and N. Now imagine that the population goes back up from its
#    current low? We would predict an increase in kappa, correct?
#    Would we still predict the same increase if the number of hunters stays low?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



