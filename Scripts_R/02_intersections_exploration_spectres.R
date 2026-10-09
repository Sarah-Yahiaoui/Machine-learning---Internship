# ============================================================
# 02 - Intersections AAOmega/SDSS et exploration des spectres
# ============================================================
# Objectif :
#   - charger les spectres AAOmega ;
#   - calculer les intersections avec l'échantillon T18 et SDSS ;
#   - visualiser quelques spectres et comparer des galaxies.
#
# Fichiers d'entrée attendus :
#   GAMA_PWuWoa.csv
#   GAMA_2556.csv
#   fichiers FITS AAOmega et SDSS cités dans le script.
#
# Ce script utilise les objets créés dans 01_construction_echantillon.R.
# ============================================================

library(dplyr)
library(FITSio)

################  Lecture des spectres avec le survey AAOmega ####################################################

Spectre_AAO <- read.csv("GAMA_PWuWoa.csv")
names(Spectre_AAO)
dim(Spectre_AAO) 

 

# Maintenant je vais étudier les intersections avec les autres survey: SDSS et l'échantillon GAMMA (7338)

# 1.Intersection avec GAMMA 7338 

save(Spectre_AAO,file = "AAOmega.RData")

length(intersect(GAMMA$CATAID, Spectre_AAO$CATAID)) 

GAMMAetAAO <- merge(GAMMA, Spectre_AAO, by='CATAID')
dim(GAMMAetAAO) 

# il y'a 3740 galxies qui appartienne à l'échantillon T18 sur 3870

# 2.Intersection avec SDSS

gamma_sdss <- read.csv("GAMA_2556.csv")
names(gamma_sdss)
dim(gamma_sdss)

GAMMA_SDSS <- merge(gamma_sdss, GAMMA, by='CATAID') 
dim(GAMMA_SDSS)

length(intersect(GAMMA_SDSS$CATAID,Spectre_AAO$CATAID))

# 3740 - 281 

SDSSetAAO <- merge(GAMMA_SDSS,Spectre_AAO, by='CATAID')
dim(SDSSetAAO)


#Taitement des spectres AAOmega, plotter un spectre et le comparer avec celui des SDSS

library(FITSio)

Sp1_AAO <- readFITS("G12_Y6_053_266.fit")
Sp1_AAO$header   # la longueur d'onde est entre  3727.88 ; 8858.87 



flux_AAO1 <- Sp1_AAO$imDat[,1]
flux_AAO1[abs(flux_AAO1) > 1e4] <- NA
flux_net1 <- !is.na(flux_AAO1)
DIM1  <- length(flux_AAO1)
DIM1

crpix_AAO1 <- Sp1_AAO$axDat$crpix[1]
crval_AAO1 <- Sp1_AAO$axDat$crval[1]
cdelt_AAO1 <- Sp1_AAO$axDat$cdelt[1]

lambda_Sp1_AAO <- crval_AAO1 + (1:DIM1 - crpix_AAO1)* 1.035935498767 
range(lambda_Sp1_AAO)  # [3816.86 ; 8769.86] 



plot(lambda_Sp1_AAO[flux_net1], flux_AAO1[flux_net1],
     type="l",
     col="black",
     main="Spectre_AAOmega(CATAID  611445)",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")
     #,ylim = c(0,200))

# Ce que je vient de plotter c'est le spectre de la galaxie sous le CATAID  611445

# Voir son taux de formation stellaire et sa masse

611445 %in% SDSSetAAO$CATAID
611445 %in% GAMMAetAAO$CATAID

which(GAMMAetAAO$CATAID == 611445)
GAMMAetAAO$CATAID[3243]
GAMMAetAAO$sSFR_0_1Gyr_best_fit[3243] # = 3.679e-10 inclu [1.770e-14,9.375e-09]

summary(GAMMAetAAO$sSFR_0_1Gyr_best_fit)

GAMMAetAAO$mass_stellar_best_fit[3243]
summary(GAMMAetAAO$mass_stellar_best_fit) # = 4959000 inclu [1.171e+06,3.134e+11]

# Je dirais plutôt un grand taux de formation stellaire et petite masse stellaire




# Je vais plotter avec la grande masse et petit taux de formation

which(GAMMAetAAO$mass_stellar_best_fit == 3.134e+11) # = 1296
GAMMAetAAO$URL[1296]
GAMMAetAAO$sSFR_0_1Gyr_best_fit[1296]

GAMMAetAAO$CATAID[1296]
Sp2_AAO <- readFITS("G15_Y4_208_085.fit")
Sp2_AAO$header
250251 %in% GAMMAetAAO$CATAID
250251 %in% SDSSetAAO$CATAID

# Cette galaxie appartient aux 2 intersections

flux_AAO2 <- Sp2_AAO$imDat[,1]
flux_AAO2[abs(flux_AAO2) > 1e4] <- NA
flux_net2 <- !is.na(flux_AAO2)
DIM2  <- length(flux_AAO2)

crpix_AAO2 <- Sp2_AAO$axDat$crpix[1]
crval_AAO2 <- Sp2_AAO$axDat$crval[1]
cdelt_AAO2 <- Sp2_AAO$axDat$cdelt[1]

lambda_Sp2_AAO <- crval_AAO2 + (1:DIM2 - crpix_AAO2)*1.035948462422

plot(lambda_Sp2_AAO[flux_net2], flux_AAO2[flux_net2],
     type="l",
     col="black",
     main="Spectre_AAOmega(CATAID 250251 )",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")

# cette galaxie a un petit taux de formation stellaire et une grande masse 

# Maintenant je vais prendre 2 galaxies dans les 2 surveys AAOmega et SDSS avec petite et grande masse et petit et grand taux de formation stellaire 
# je vais prendre une galaxie qui appartient à lintersection des 3 ech ( la ou il ya 281 galaxies)

# Je commence par AAOmega avec une petite masse et grand taux de formation

summary(SDSSetAAO$sSFR_0_1Gyr_best_fit)
which(SDSSetAAO$sSFR_0_1Gyr_best_fit == 8.441e-09)
SDSSetAAO[129,]
SDSSetAAO$sSFR_0_1Gyr_best_fit[126]

# C'est la galaxie sous  CATAID 271562 avec le grand taux de formation et de masse 18500000 plutot petite

summary(SDSSetAAO$mass_stellar_best_fit)


# Le plot avec SDSS

Sp2_sdss <- readFITS("spec-0513-51989-0054 (2).fit")

flux_Sp2_sdss <- Sp2_sdss$imDat[,1]
flux_Sp2_sdss[abs(flux_Sp2_sdss) > 1e4] <- NA
flux_Sp2_sdss_net <- !is.na(flux_Sp2_sdss)
DIM_sdss2 <- length(flux_Sp2_sdss)

crpix_Sp2_sdss <- Sp2_sdss$axDat$crpix[1]
crval_Sp2_sdss <- Sp2_sdss$axDat$crval[1]
cdelt_Sp2_sdss <- Sp2_sdss$axDat$cdelt[1]

LOGlambda_Sp2_sdss <- crval_Sp2_sdss + (1:DIM_sdss2 - crpix_Sp2_sdss)*cdelt_Sp2_sdss
lambda_Sp2_sdss<- 10^LOGlambda_Sp2_sdss
lambda_Sp2_sdss[2] - lambda_Sp2_sdss[1]
par(mfrow=c(1,2), oma=c(0,0,2,0))

plot(lambda_Sp2_sdss[flux_Sp2_sdss_net], flux_Sp2_sdss[flux_Sp2_sdss_net],
     type="l",
     col="blue",
     main="Spectre_SDSS",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")

# Le plot avec AAOmega

Sp3_AAO <- readFITS("G12_Y6_070_258.fit")
Sp3_AAO$header
Sp3_AAO$hdr
flux_AAO3 <- Sp3_AAO$imDat[,1]
flux_AAO3[abs(flux_AAO3) > 1e4] <- NA
flux_net3 <- !is.na(flux_AAO3)
DIM3  <- length(flux_AAO3)


crpix_AAO3 <- Sp3_AAO$axDat$crpix[1]
crval_AAO3 <- Sp3_AAO$axDat$crval[1]
cdelt_AAO3 <- Sp3_AAO$axDat$cdelt[1]

lambda_Sp3_AAO <- crval_AAO3 + (1:DIM3 - crpix_AAO3)* 1.035935498767
dt <- lambda_Sp3_AAO[2] -lambda_Sp3_AAO[1]
dt
plot(lambda_Sp3_AAO[flux_net3], flux_AAO3[flux_net3],
     type="l",
     col="black",
     main="Spectre_AAOmega(CATAID)",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")
mtext("Galaxie 271562 avec une grand taux de formation", outer = TRUE, line = -1)

# plot d'une galaxie avec une grande masse

summary(SDSSetAAO$mass_stellar_best_fit)
which(SDSSetAAO$mass_stellar_best_fit == 3.134e+11) #112
SDSSetAAO[112,] 

# La galaxie est sous le cataid 250251

# Plot sous SDSS

Sp3_sdss <- readFITS("spec-4030-55634-0118.fit")

flux_Sp3_sdss <- Sp3_sdss$imDat[,1]
flux_Sp3_sdss[abs(flux_Sp3_sdss) > 1e4] <- NA
flux_Sp3_sdss_net <- !is.na(flux_Sp3_sdss)
DIM_sdss3 <- length(flux_Sp3_sdss)

crpix_Sp3_sdss <- Sp3_sdss$axDat$crpix[1]
crval_Sp3_sdss <- Sp3_sdss$axDat$crval[1]
cdelt_Sp3_sdss <- Sp3_sdss$axDat$cdelt[1]

LOGlambda_Sp3_sdss <- crval_Sp3_sdss + (1:DIM_sdss3 - crpix_Sp3_sdss)*cdelt_Sp3_sdss
lambda_Sp3_sdss<- 10^LOGlambda_Sp3_sdss

par(mfrow=c(1,2), oma=c(0,0,2,0))

plot(lambda_Sp3_sdss[flux_Sp3_sdss_net], flux_Sp3_sdss[flux_Sp3_sdss_net],
     type="l",
     col="blue",
     main="Spectre_SDSS",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")

# Le plot avec AAOmega

Sp4_AAO <- readFITS("G15_Y4_208_085.fit")
Sp4_AAO$header

flux_AAO4 <- Sp4_AAO$imDat[,1]
flux_AAO4[abs(flux_AAO4) > 1e4] <- NA
flux_net4 <- !is.na(flux_AAO4)
DIM4  <- length(flux_AAO4)


crpix_AAO4 <- Sp4_AAO$axDat$crpix[1]
crval_AAO4 <- Sp4_AAO$axDat$crval[1]
cdelt_AAO4 <- Sp4_AAO$axDat$cdelt[1]

lambda_Sp4_AAO <- crval_AAO4 + (1:DIM4 - crpix_AAO4)*1.035948


##########  Test sur les valeurs de la longueur d'onde #############

range(lambda_Sp4_AAO) # =  [ 3815.534 ; 8768.534 ]
Sp4_AAO$header        # = [ WMIN ; WMAX ] = [ 3726.53 ; 8857.58 ]



plot(lambda_Sp4_AAO[flux_net4], flux_AAO4[flux_net4],
     type="l",
     col="black",
     main="Spectre_AAOmega",
     xlab="Longueur d'onde (Å)",
     ylab="Flux")

mtext("Galaxie 250251 avec une grande masse", outer = TRUE, line = -1)
250251 %in% GAMMAetAAO$CATAID
which(GAMMAetAAO$CATAID==250251)
GAMMAetAAO$Z_TONRY[1296]

# Observation : les spectres plottés sont bon et se ressemble à ce qu'il ya dans la base du site GAMA les png
