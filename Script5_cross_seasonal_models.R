# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script runs a 'cross-seasonal' model to estimate summer and winter
# survival rates of adult female mallards marked in the North American midcontinent
# 
# ~ Thomas Riecke 
# thomasvanceriecke@gmail.com
# thomas.riecke@umontana.edu
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x <- url("https://github.com/thomasriecke/DUCKS9/blob/main/mallard_data.RData?raw=true")
load(x)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-seasonal model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(jagsUI)
# Specify model in BUGS language
sink("cross_seasonal.jags")
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
    
    sigma.phi ~ dgamma(1,1)
    sigma.omega ~ dgamma(1,1)
    sigma.kappa ~ dgamma(1,1)
    
    tau.phi = 1/(sigma.phi * sigma.phi)
    tau.omega = 1/(sigma.omega * sigma.omega)
    tau.kappa = 1/(sigma.kappa * sigma.kappa)

    for(t in 1:n.years){
    
      epsilon.phi[t] ~ dnorm(0, tau.phi)
      logit(phi[t]) = alpha.phi + epsilon.phi[t] 
      
      epsilon.kappa[t] ~ dnorm(0, tau.kappa)
      logit(kappa[t]) = alpha.kappa + epsilon.kappa[t]
      
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
                  rho.alpha = rho.alpha, rho.beta = rho.beta)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f','h','phi','omega','rho','kappa','c',
                'sigma.h','sigma.omega','sigma.phi',
                'epsilon.h','epsilon.omega','epsilon.phi')

nc <- 4
nt <- 25
ni <- 10000
nb <- 2500


# Call JAGS from R (BRT )
Sys.time()
m <- jags(jags.data, inits, parameters, "cross_seasonal.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A few disclaimers:
#
# We didn't constrain survivals and h fully effectively, 
# i.e., in this model, it's possible for harvest to be 0.1 and winter survival to be 0.91!
# for further reading about how to control for this type of issue see:
#
# Ergon et al. (2018) Methods in Ecology and Evolution
# Nater et al. (2020) Journal of Animal Ecology
# Riecke et al. (2022a, 2022b) Journal of Animal Ecology
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
greens <- colorRampPalette(c('lightgreen','forestgreen'))

par(mar = c(5.1,5.1,2.1,2.1), family = 'serif')
library(vioplot)
vioplot(m$sims.list$phi, names = seq(start,end), col = greens(n.years),
        las = 1, ylab = expression(Summer~survival~(phi)), cex.lab = 2,
        xlab = 'Year', drawRect = F)

vioplot(m$sims.list$omega, names = seq(start,end-1), col = greens(n.years),
        las = 1, ylab = expression(Winter~survival~without~harvest~(omega)), cex.lab = 2,
        xlab = 'Year', drawRect = F)

vioplot(m$sims.list$kappa, names = seq(start,end), col = greens(n.years),
        las = 1, ylab = expression(Hunting~mortality~probability~(h)), cex.lab = 2,
        xlab = 'Year', drawRect = F)

vioplot(m$sims.list$h, names = seq(start,end), col = greens(n.years),
        las = 1, ylab = expression(Harvest~probability~(h)), cex.lab = 2,
        xlab = 'Year', drawRect = F)

vioplot(m$sims.list$rho, names = seq(start,end), col = greens(n.years),
        las = 1, ylab = expression(Band~reporting~probability~(rho)), cex.lab = 2,
        xlab = 'Year', drawRect = F)

vioplot(m$sims.list$S, names = seq(start,end-1), col = greens(n.years),
        las = 1, ylab = expression(Survival~(S)), cex.lab = 2,
        xlab = 'Year', drawRect = F)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions:
# 
# 1a) Is annual survival (S) more closely associated with summer survival (phi)
#    or hunting mortality (kappy)? Why?
#
plot(m$q50$S ~ m$q50$phi[2:n.years], pch = 21, bg = greens(n.years-1),
     ylab = 'Annual survival', xlab = 'Summer survival', cex = 1.5, cex.lab = 1.5)
plot(m$q50$S ~ m$q50$kappa[1:(n.years-1)], pch = 21, bg = greens(n.years-1),
     ylab = 'Annual survival', xlab = 'Hunting mortality', cex = 1.5, cex.lab = 1.5)
# 
# 1b) Why? See Hoekman et al. (2002; JWM) and Koons et al. (2016; Ecology Letters)
#
# 2) What happens to our estimates if we use a 'weaker' prior for crippling loss probability?
#    What might that prior look like? A stronger prior?
#
# 3) How might we incorporate a covariate (e.g., Ponds) to help estimate summer survival?
# 
# A key take-home here is that these types of models (cross-seasonal models) are incredibly 
# useful. However, as the complexity of our models increases, they require additional
# assumptions. Post-season banding is important for estimating these parameters.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




