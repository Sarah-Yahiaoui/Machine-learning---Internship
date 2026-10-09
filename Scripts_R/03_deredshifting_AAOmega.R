# ============================================================
# 03 - Préparation des spectres AAOmega et dé-redshifting
# ============================================================
# Objectif :
#   - lire les fichiers FITS AAOmega ;
#   - reconstruire les longueurs d'onde observées ;
#   - corriger les longueurs d'onde du redshift ;
#   - déterminer une plage commune de longueurs d'onde.
#
# Ce script utilise les objets créés dans les scripts 01 et 02.
# ============================================================

library(FITSio)

################# ETAPE 2: DE-REDSHIFTER ########################

# Je commence par déredshifter pour ensuite plotter et voir les valeurs min et max de la longueur d'onde comme ce qu'il ya dans SDSS 

dim(Spectre_AAO)
names(Spectre_AAO)

# Pour dé-redshifer on va utiliser que ce qui appartient à l'échantillon 
# c'est à dire les 3740 (avec 281 qui appartienne à sdss aussi)

dim(GAMMAetAAO)
names(GAMMAetAAO)


files_AAOmegaetGAMMA <- list.files("FITS_AAOmegaetGAMMA(3740)",pattern = "\\.fit$", full.names = TRUE)
length(files_AAOmegaetGAMMA)  # = 3737



# 1. je vais créer une liste des spectres, pour ensuite les corrigés dans une autre liste 
# je vais créer une boucle qui va parcourir tout le fichier 3740 (les spectres) et leurs calculer leurs longueurs d'onde 

lambdas_obs <- list()  # création des listes vides 
flux_AAO <- list()

# la boucle

for (i in 1:length(files_AAOmegaetGAMMA)){
  
  # je doit d'abord lire les fochiers en  fits
  
  spectres_AAO <- readFITS(files_AAOmegaetGAMMA[i])
  
  # je vais extraire de chaque spectres les cdelt crval et crpix pour calculer lambda
  
  crval <- spectres_AAO$axDat$crval[1]
  crpix <- spectres_AAO$axDat$crpix[1]
  cdelt_1 <- grep("CD1_1", spectres_AAO$header, value = TRUE)
  cdelt <- as.numeric(strsplit(strsplit(cdelt_1, "=")[[1]][2], "/")[[1]][1])
  
  # Extraire le flux de chaque spectre
  
  spectre_AAO_flux <- spectres_AAO$imDat[,1]
  spectre_AAO_flux[abs(spectre_AAO_flux) >1e4] <- NA
  N_flux <- length(spectre_AAO_flux)
  
  # Calcul de la longueur d'onde observée
  
  lambda_AAO_obs <- crval +(1:N_flux - crpix)*cdelt
  
  # je vais faire la sauvegarde de tout les lambda observé dans la liste déjà crée
  
  lambdas_obs[[i]] <- lambda_AAO_obs
  flux_AAO[[i]] <- spectre_AAO_flux
  cat(i, "/", length(files_AAOmegaetGAMMA), " traité\n")  # juste pour les afficher comment il fait le traitement 
}


# 1.ere Etape : Correction des longueur d'onde 

length(files_AAOmegaetGAMMA)

lambdas_corr_AAO <- list()

# test sur le redshift du spectre

sp <- readFITS(files_AAOmegaetGAMMA[1])
sp$header 
Z_txt <- sp$header[grep("^Z\\s*=", sp$header)]

Z <- as.numeric(
  sub(".*=\\s*([0-9.]+).*", "\\1", Z_txt)
)

Z # 0.02932

Z_AAO <- numeric(length(files_AAOmegaetGAMMA))

# Boucle sur tout les spectres

for(i in seq_along(lambdas_obs)) {
  
  sp <- readFITS(files_AAOmegaetGAMMA[i])
  
  Z_txt <- sp$header[grep("^Z\\s*=", sp$header)]
  
  Z <- as.numeric(sub(".*=\\s*([0-9.]+).*", "\\1", Z_txt)
  )
  
  Z_AAO[i] <- Z
  
  lambda_corr <- lambdas_obs[[i]] / (1 + Z)
  
  lambdas_corr_AAO[[i]] <- lambda_corr # c'est une list [[i]]
}



# la longueur d'onde minimale est la valeur max des minima de tous les spectres
# la longueur d'onde maximale est la valeur min des maxima e tous les spectres 
# le pas c'est le min de tous les pas de tous les spectres 

range(sapply(lambdas_corr_AAO, max))

mins <- sapply(lambdas_corr_AAO,min)

lambda_min <- max(mins)
lambda_min 


Maxs <- sapply(lambdas_corr_AAO,max)

lambda_max <- min(Maxs)
lambda_max   

# la fonction diff() prend la différence entre chanque valeur consécutives et le min pour prendre le min entre elle


steps <- sapply(lambdas_corr_AAO, function(x) min(diff(x)))

step_min <- min(steps)
step_min 

#===================================================
# Proportion des spectres 
#=====================================================

lambda_max_all <- sapply(lambdas_corr_AAO, max)
n_6000 <- sum(lambda_max_all < 8350)
n_6000 # 1   

# Effectivement il n'y a qu'un seul spectre qui termine a 6230
# identification de ce spectre

GAMMAetAAO$URL[which(lambda_max_all <= 6230)]
GAMMAetAAO$CATAID[which(lambda_max_all <= 6230)] 

n_7000 <- sum(lambda_max_all < 7000)
n_7000
n_8835 <- sum(lambda_max_all <= 8834 & lambda_max_all >= 8500 )
n_8835

n_6000_7000 <- sum(lambda_max_all < 7000 & lambda_max_all> 6230)
n_6000_7000 # 2

n_7000_8000 <- sum(lambda_max_all < 8000 & lambda_max_all > 7000)
n_7000_8000

n_8000 <- sum(lambda_max_all > 8350)
n_8000 # 3

n_9000 <- sum(lambda_max_all < 8500 & lambda_max_all >=8500)
n_9000


range(lambdas_corr_AAO)

# il y'en un spectre dont la longueur d'onde maximale fini à 6230 (le max de la grille)
# il ya 2 spectre entre 6230 et 7000 et le reste leur max fini tous après 8000A

i_bad <- which.min(lambda_max_all) # 3593 c'est le numero du spectre qui donne la valeur minimale de la grille 
GAMMAetAAO$CATAID[i_bad] 

idx_remove <- which(lambda_max_all < 8350)

idx_remove
data.frame(
  indice = idx_remove,
  CATAID = GAMMAetAAO$CATAID[idx_remove],
  lambda_max = lambda_max_all[idx_remove]
)

# j vais enlever les 3 spectres qui sont trop petit :longuer d'onde max petite que les autres ce qui fait qu'elle ecrase de l'information


CATAID_recherche <- c(485777, 543535, 3852250)

files <- list.files("fits_AAOMEGA",
                    pattern="\\.fit$",
                    full.names=TRUE)

for(id in CATAID_recherche){
  
  trouve <- FALSE
  
  for(f in files){
    
    sp <- readFITS(f)
    
    ligne <- grep("CATAID", sp$header, value=TRUE)
    
    cataid_fit <- as.numeric(
      strsplit(
        strsplit(ligne, "=")[[1]][2],
        "/"
      )[[1]][1]
    )
    
    if(cataid_fit == id){
      
      cat("\nCATAID :", id, "\n")
      cat("Fichier :", basename(f), "\n")
      
      trouve <- TRUE
      break
    }
  }
  
  if(!trouve){
    cat("\nCATAID", id, "non trouvé.\n")
  }
}
