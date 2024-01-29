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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-seasonal model with Lincoln estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(jagsUI)
# Specify model in BUGS language
sink("cross_seasonal_Lincoln.jags")
cat("
    model {


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Lincoln estimates
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (t in 1:n.years){
    
      N[t,1] = N[t,2]/phi[t]
      N[t,2] = H[t]/h[t]
      
      
      # The hashed out approach would apply additional uncertainty to
      # the estimates of H (Harvest)
      # N[t,1] = N[t,2]/phi[t]
      # N[t,2] = Harvest[t]/h[t]
      # Harvest[t] ~ dnorm(H[t], tauH[t])
      # sigmaH[t] = H[t] * 0.15
      # tauH[t] = 1/(sigmaH[t] * sigmaH[t])
      
    }


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
                  rho.alpha = rho.alpha, rho.beta = rho.beta, H = H)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('S','f','h','phi','omega','rho','kappa','c',
                'sigma.h','sigma.omega','sigma.phi',
                'epsilon.h','epsilon.omega','epsilon.phi','N')

nc <- 4
nt <- 25
ni <- 10000
nb <- 2500


# Call JAGS from R (BRT )
Sys.time()
m6 <- jags(jags.data, inits, parameters, "cross_seasonal_Lincoln.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()




greens <- colorRampPalette(c('lightgreen','darkgreen'))
plot(m6$q50$N[,1] ~ seq(start,end), ylim = c(2000000,12000000),
     ylab = 'Abundance (N)', xlab = 'Year', cex.lab = 2)
arrows(seq(start,end), m6$q2.5$N[,1], seq(start,end), m6$q97.5$N[,1],
       length = 0, lty = 1)
points(m6$q50$N[,1] ~ seq(start,end), cex = 2, type = 'b',
       pch = 21, bg = greens(n.years))

arrows(seq(start+0.25,end+0.25), m6$q2.5$N[,2], seq(start+0.25,end+0.25), m6$q97.5$N[,2],
       length = 0, lty = 1)
points(m6$q50$N[,2] ~ seq(start+0.25,end+0.25), cex = 2, type = 'b',
       pch = 22, bg = greens(n.years))

reds <- colorRampPalette(c('orange','red'))
points(c(N*1000000) ~ seq(start,end), bg = reds(n.years), 
       type = 'b', pch = 21, cex = 1.5)

legend('topright', legend = c('Lincoln spring', 'Lincoln fall', 'WBPHS'),
       pch = c(21,22,21), cex = 1.25,
       pt.bg = c('darkgreen','darkgreen','red'), bty = 'n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Questions
# 1) What happens to the confidence intervals of Lincoln estimates 
#    if you apply additional uncertainty to the estimates of H (harvest)?
#    see the hashed out code in the model
#
# 2) What happens if band-reporting probabilities are different than our 
#    assumptions? Try changing the prior for band-reporting probabilities
#    to ~0.3 by adding the following code and re-running the model
#    rho.alpha <- rep(30, n.years)
#    rho.beta <- rep(70, n.years)
#
#    What if we were more certain of this new estimate? Try this code:
#    rho.alpha <- rep(90, n.years)
#    rho.beta <- rep(210, n.years)
#    and re-run the model
# 
# 3) Reconstruct the plot above minus the fall estimates (remove N[,2]).
#    Do Lincoln estimates and WBPHS estimates generally agree?
#    Given the assumptions of Lincoln estimates (see Alisauskas papers)
#    what might drive divergence in these two estimates of abundance
#
# 4) Why are spring estimates of abundance greater than fall estimates?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



