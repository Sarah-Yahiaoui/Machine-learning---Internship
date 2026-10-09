# ============================================================
# 01 - Construction de l'échantillon T18 et paramètres physiques
# ============================================================
# Objectif :
#   - charger la base GAMA ;
#   - appliquer les filtres de l'échantillon ;
#   - joindre les distances ;
#   - calculer les paramètres physiques et les standardiser.
#
# Fichiers d'entrée attendus dans le répertoire de travail :
#   GAMA_DATA_7556.csv
#   DistancesFramesv12.fits
#
# Remarque : vérifier les effectifs après chaque filtre.
# ============================================================

library(dplyr)
library(FITSio)

#########################################################################################################"""
# Dans ce script, il y'a l'échantillon de T18 avec les 5 paramètres physiques et les spectres AAOmega
################################################################################################"


gamma <- read.csv("GAMA_DATA_7556.csv") 

# Maintenant je doit les filtrer pour passer à 7338

# 1. supprimer les valeurs manquantes

nrow(gamma)
names(gamma)
library(dplyr)
df1 <- na.omit(gamma)
nrow(df1)

#2. Passage de 7540 galaxie a 7516 
library(dplyr)
df2 <- df1 %>%
  filter(!HUBBLE_TYPE %in% c("Star", "Artifact"))
nrow(df2)
names(df2)


#3. Applications des conditions aux limites tels que T18



# c'est à partir de là qu'on fait la conversion

nrow(df2) 


library(FITSio)
DM <- readFITS("DistancesFramesv12.fits")
DM$colNames ### Pour voir les colonne CATAID ET DM pour les joindre a la table apres 
DM12 <- as.data.frame(DM$col)
colnames(DM12) <- DM$colNames
head(DM12)
names(DM12)
nrow(DM12)

gamma_DM <- left_join(df2, DM12, by='CATAID') 
sum(is.na(gamma_DM$DM_70_30_70))

names(gamma_DM)
nrow(gamma_DM)

## maintenant je vais convertir en Kpc

library(dplyr)

gamma_DM <- gamma_DM %>%
  mutate(
    DL = 10^((DM_70_30_70 - 25)/5),
    DA = DL / (1 + Z_TONRY.x)^2
  )

names(gamma_DM)
nrow(gamma_DM)

summary(gamma_DM$Z_TONRY.x) 
summary(gamma_DM $Z_TONRY.y)

# ici Z_TONRY.x=Z_TONRY.y

gamma_DM <- gamma_DM %>%
  mutate(
    Rhalf_kpc = GALRE_r * DA * 1000 / 206265
  )
gamma_DM <- gamma_DM %>%
  filter(
    mass_stellar_best_fit > 0,
    GALINDEX_r > 0,
    sSFR_0_1Gyr_best_fit > 0
  )

gamma_DM <- gamma_DM %>%
  mutate(
    logRhalf = log10(Rhalf_kpc),
    logMstar = log10(mass_stellar_best_fit),
    logn = log10(GALINDEX_r),
    logsSFR = log10(sSFR_0_1Gyr_best_fit)
  )
nrow(gamma_DM)

GM <- gamma_DM %>%
  filter(
    logMstar >= 6 & logMstar <= 12,
    uminusr >= 0.3 & uminusr <= 2.7,
    logn >= (-0.6) & logn <= 1.2,     #log10(GALRE_r)
    logRhalf >= (-1.0) & logRhalf <= 1.5,   # R_1SUR2
    logsSFR >= (-14) & logsSFR <=  (-8)     
  )



# Ce que je doit faire c'est d'enlever les variables qui était dans la base DMV12 pour avoir la même dimension
GAMMA  <- GM[ ,
              !names(GM) %in% c(
                "RA","DEC","DM_70_30_70","DM_100_30_70","DM_70_25_75","DM_100_25_75","Q","Z_HELIO","Z_LG",  
                "Z_CMB","Z_TONRY.x" , "Z_TONRY_HIGH","Z_TONRY_LOW","DL","DA"   
              )
]
names(GAMMA)[names(GAMMA) == "Z_TONRY.y"] <- "Z_TONRY"


## Observation : la remarque est juste, je doit normaliser toute mes valeurs en log

GAMMA$logRhalf   <- scale(GAMMA$logRhalf)
GAMMA$logn      <- scale(GAMMA$logn)
GAMMA$logMstar <- scale(GAMMA$logMstar)
GAMMA$logsSFR    <- scale(GAMMA$logsSFR)
uminusr_Z <- scale(GAMMA$uminusr)
GAMMA <-cbind(GAMMA,uminusr_Z )


#################################################################################################################################################################

# ce qui est en log dans GAMMA est _z pour z-score pour GAMA_DATA_7338
