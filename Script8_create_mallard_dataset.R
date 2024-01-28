# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this script reads and formats mallard banding, encounter, and
# harvest data to estimate survival, natural mortality, and hunting mortality
# probabilities, crippling loss, band reporting probability, and abundance
# for each age (HY/AHY) and sex (M/F) class
# 
# ~ Thomas Riecke 
# thomasvanceriecke@gmail.com
# thomas.riecke@umontana.edu
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(plyr)
library(readxl)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read and format banding data
# note you'll have to change the directory to wherever you placed the 
# file 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read.csv('E:/O_Final/Ducks/banding_data/2022/NABBP_2022_grp_08.csv', header = TRUE)
table(dat$EVENT_TYPE)
rel <- subset(dat, EVENT_TYPE == 'B')
enc <- subset(dat, EVENT_TYPE == 'E')
rm(dat)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of years
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <- 1974
end <- 1994
n.years <- end - start + 1
n.seasons <- 2



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# important banding/encounter codes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(rel$BIRD_STATUS)
table(rel$ISO_SUBDIVISION)
table(rel$EXTRA_INFO)
table(rel$BAND_TYPE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# acceptable band types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://www.sciencebase.gov/catalog/item/632b2d7bd34e71c6d67bc161
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ok.bands <- c('01','04','08','11','12','13','18','21','41','51','53','81','82','W1')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# acceptable 'extra info'
# 00 = federal numbered metal band only
# 04 = control band (for use in conjunction with reward band studies only)
# 70 = captured by spotlighting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ok.info <- c('00', '04', '70')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# acceptable release locations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://www.sciencebase.gov/catalog/item/632b2d7bd34e71c6d67bc161
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ok.release.locations <- c('CA-AB','CA-SK','CA-MB',
                          'US-ND','US-SD','US-MT','US-MN',
                          'US-NE','US-KS','US-OK','US-TX',
                          'US-IA','US-MO','US-AR','US-LA')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format pre-season banding releases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(rel)
rel <- subset(rel, ISO_SUBDIVISION %in% ok.release.locations)
rel <- subset(rel, EVENT_YEAR >= start & EVENT_YEAR <= end)
rel <- subset(rel, BIRD_STATUS == '3')
rel <- subset(rel, BAND_TYPE %in% ok.bands)
rel <- subset(rel, EXTRA_INFO %in% ok.info)
table(rel$EVENT_YEAR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format pre-season banding encounters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(enc$HOW_OBTAINED)
table(enc$WHO_OBTAINED)
table(enc$BAND_TYPE)
table(enc$ISO_SUBDIVISION)

enc <- subset(enc, EVENT_YEAR >= start & EVENT_YEAR <= (end + 1))
enc <- subset(enc, HOW_OBTAINED == '1')
enc <- subset(enc, WHO_OBTAINED == '21')
enc <- subset(enc, BAND_TYPE %in% ok.bands)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# handle recovery months (i.e., make Jan 1962 recovery correspond to 1961 release)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(enc$EVENT_MONTH)
table(enc$EVENT_YEAR)
ok.months <- c(1,2,9,10,11,12,92,93,94)
enc <- subset(enc, EVENT_MONTH %in% ok.months)
enc$EVENT_YEAR[which(enc$EVENT_MONTH %in% c(1,2,92))] <- enc$EVENT_YEAR[which(enc$EVENT_MONTH %in% c(1,2,92))] - 1
table(enc$EVENT_YEAR)
enc <- subset(enc, EVENT_YEAR >= start & EVENT_YEAR <= end)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# join encounters to releases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enc <- subset(enc, select = c('BAND','EVENT_YEAR'))
names(enc) <- c('BAND','R_YEAR')
dat <- join(rel, enc, type = 'left')

table(dat$R_YEAR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# eliminate impossible encounters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$ENC_TIME <- dat$R_YEAR - dat$EVENT_YEAR
table(dat$ENC_TIME)
dat <- subset(dat, ENC_TIME >= 0 | is.na(ENC_TIME))
dat <- subset(dat, ENC_TIME <= 25 | is.na(ENC_TIME))
table(dat$ENC_TIME)

rm(enc, rel, ok.bands, ok.months, ok.release.locations, ok.info)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format m-arrays
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(dat$AGE_CODE)
table(dat$SEX_CODE)

rel <- matrix(0, n.years, n.seasons)
marr <- array(0, dim = c(n.years, n.years + 1, n.seasons))

dat.af <- subset(dat, AGE_CODE %in% c(1,5,6,8) & SEX_CODE == '5')

pre.months <- c(7,8,9)
post.months <- c(1,2,3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# post-season adult females: m[,,1], rel[,1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- subset(dat.af, EVENT_MONTH %in% post.months)

for (i in 1:nrow(tmp)){
  rel[tmp$EVENT_YEAR[i] - start + 1, 1] <- rel[tmp$EVENT_YEAR[i] - start + 1, 1] + 1
  
  if(!is.na(tmp$R_YEAR[i])){
    marr[tmp$EVENT_YEAR[i] - (start - 1), tmp$R_YEAR[i] - (start - 1), 1] <- 
      marr[tmp$EVENT_YEAR[i] - (start - 1), tmp$R_YEAR[i] - (start - 1), 1] + 1
  }
}

marr[,n.years+1,1] <- rel[,1] - rowSums(marr[,1:n.years,1])



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pre-season adult females: m[,,2], rel[,2]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- subset(dat.af, EVENT_MONTH %in% pre.months)

for (i in 1:nrow(tmp)){
  rel[tmp$EVENT_YEAR[i] - start + 1, 2] <- rel[tmp$EVENT_YEAR[i] - start + 1, 2] + 1
  
  if(!is.na(tmp$R_YEAR[i])){
    marr[tmp$EVENT_YEAR[i] - (start - 1), tmp$R_YEAR[i] - (start - 1), 2] <- 
      marr[tmp$EVENT_YEAR[i] - (start - 1), tmp$R_YEAR[i] - (start - 1), 2] + 1
  }
}

marr[,n.years+1,2] <- rel[,2] - rowSums(marr[,1:n.years,2])


# View m-arrays
marr[,,1]
marr[,,2]


par(mfrow = c(1,2))
barplot(rel[,1], ylim = c(0,20000), ylab = 'Post-season releases')
barplot(rel[,2], ylim = c(0,20000), ylab = 'Post-season releases')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# clean up
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(dat, post.months, pre.months, dat.af,tmp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# additional data...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read_xlsx('E:/O_Final/Duck_Symposium_9/band_recovery_workshop/data.xlsx')
names(dat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Abundance and pond data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N <- dat$BPOP[(start-1960):(end-1960)]/1000
P <- as.numeric(dat$Ponds[(start-1960):(end-1960)])/1000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Harvest data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
H <- dat$Htotal[(start-1960):(end-1960)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Duck stamp data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D <- dat$Stamps[(start-1960):(end-1960)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reporting rates (moment matching)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rho.mu <- rep(0.38, n.years)
rho.sd <- rep(0.02, n.years)

rho.alpha = ((1-rho.mu)/(rho.sd*rho.sd) - (1/rho.mu)) * rho.mu^2
rho.beta = rho.alpha * (1/rho.mu - 1)
rm(dat)

save.image("E:/O_Final/Duck_Symposium_9/band_recovery_workshop/final_scripts/mallard_data.RData")

load("E:/O_Final/Duck_Symposium_9/band_recovery_workshop/final_scripts/mallard_data.RData")





