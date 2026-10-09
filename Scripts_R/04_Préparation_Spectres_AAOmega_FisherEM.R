#====================================================================================
# AAOmega rééchantillonnées sur leur intervalle 
#================================================================================

# Ici les longueures d'onde sont corrigés 

length(lambdas_corr_AAO)
range(lambdas_corr_AAO)  # = 3518.262 8834.345


#================================================================
# Le resample avec la fonction approx
#================================================================

# Définir la grille de longueur d'onde 

lambda_AAOgrid <- seq(3724.337, 8350, by=0.7283717)

flux <- flux_AAO[[i]]

ok <- is.finite(flux)

lambda <- lambdas_corr_AAO[[i]]

flux_sans_NA <- flux[ok]

lambda_sans_NA <- lambda[ok]

spectre_resampled_AAO_approx <- function(lambda_sans_NA, flux_sans_NA, lambda_AAOgrid){
  approx(x=lambda_sans_NA,
         y=flux_sans_NA,
         xout=lambda_AAOgrid,
         rule=1)$y}

# Premier plot pour voir si ça colle bien 

flux_AAO_resamp <- spectre_resampled_AAO_approx(lambda_sans_NA, flux_sans_NA, lambda_AAOgrid)

pdf("Spectre_resample2.pdf", width = 18, height = 7)

# 1 ligne, 2 colonnes
par(mfrow = c(1, 2),
    mar = c(5, 5, 3, 1))   # marges

# Spectre complet
plot(lambda_sans_NA, flux_sans_NA,
     type = "l",
     col = "red",
     lwd = 1,
     xlab = "lambda Å",
     ylab = "Flux",
     main = "(a) Spectre complet")
   
lines(lambda_AAOgrid, flux_AAO_resamp,
      col = "blue",
      lwd = 1
      )

legend("topright",
       legend = c("Original", "Resamplé"),
       col = c("red", "blue"),
       lwd = 2,
       bty = "n")

# Zoom sur Hα
plot(lambda_sans_NA, flux_sans_NA,
     type = "l",
     xlim = c(6400, 6600),
     col = "red",
     lwd = 1,
     xlab = "lambda Å",
     ylab = "Flux",
     main="(b) Zoom entre 6400 et 6600 Å")

lines(lambda_AAOgrid, flux_AAO_resamp,
      col = "blue",
      lwd = 1
      )

legend("topright",
       legend = c("Original", "Resamplé"),
       col = c("red", "blue"),
       lwd = 2,
       bty = "n")

dev.off()

# Resample tout les spectres 

flux_AAO_resampled_approx <- vector("list",length(flux_AAO)) 


for (i in seq_along(flux_AAO)){
  lambda <- lambdas_corr_AAO[[i]]
  flux <- flux_AAO[[i]]
  
  ok <- is.finite(flux)
  flux_sans_NA <- flux[ok]
  lambda_sans_NA <- lambda[ok]
  
  flux_AAO_resampled_approx[[i]] <- spectre_resampled_AAO_approx(lambda_sans_NA, flux_sans_NA,lambda_AAOgrid)
}

#========================================================

# Normalisation des flux

#==============================================================

intervalle_norma_AAO <- lambda_AAOgrid <= 5800 & lambda_AAOgrid >= 5600



flux_AAOmoyen <- mean(flux_AAO_resampled_approx[[i]][intervalle_norma_AAO], na.rm = TRUE)

flux_AAOnorm <- flux_AAO_resampled_approx[[i]] / flux_AAOmoyen

summary(flux_AAOnorm)  # la moyennne est autour de 1.0167


flux_AAOnorm <- lapply(
  flux_AAO_resampled_approx,
  function(f){
    flux_AAOmoyen <- mean(f[intervalle_norma_AAO], na.rm = TRUE)
    
    f / flux_AAOmoyen
  }
)

range(which(lambda_AAOgrid >= 5600 & lambda_AAOgrid <= 5800)) # l'intervalle_norm =[1397 1668]


dir.create("Spectres_AAO_norm")

for (i in seq_along(flux_AAOnorm)){
  
  spectre <- data.frame(
    
    lambda=lambda_AAOgrid,
    
    flux=flux_AAOnorm[[i]]
  )
  
  cat_id <- GAMMAetAAO$CATAID[i]
  
  
  file_name  <- paste0("Spectres_AAO_norm/spectre_",cat_id,".txt")
  
  # le cataid dans la premiere ligne du spectre
  writeLines(as.character(cat_id), file_name)
  
  write.table(
    spectre,
    file = file_name,
    sep="\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    append = TRUE
  )
}

#=================================================================
#   Création de la matrice de spectres 
#====================================================================

directory <- "Spectres_AAO_norm"

files <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)

length(files)


# Lecture du premier fichier pour connaitre le nombre de longueures d'ondes

tmp <- read.table(files[1], skip=1, header=FALSE)

n_lambda <- nrow(tmp) # ca doit etre 6351 
n_lambda # 

AAOGrid_CATAID <- matrix(nrow=length(files), ncol=1)
AAOGrid_spec   <- matrix(NA,
                             nrow=length(files),
                             ncol=n_lambda)

for(i in seq_along(files)){
  
  fileName <- files[i]
  # CATAID = première ligne
  cat_id <- readLines(fileName, n=1)
  
  # Spectre
  dat <- read.table(fileName,
                    skip=1,
                    header=FALSE)
  
  
  lambda <- as.numeric(dat[,1])
  flux    <- as.numeric(dat[,2])
  if(i==1848){
    print(length(flux))
    print(n_lambda)
    print(sum(is.na(flux)))
    print(which(is.na(flux)))
  }
  
  AAOGrid_spec[i,] <- flux
  AAOGrid_CATAID[i, 1] <- cat_id
}


AAOGrid_lbd <- lambda

AAOGrid_CATAID <- gsub("\"", "", AAOGrid_CATAID)

rownames(AAOGrid_spec) <- AAOGrid_CATAID

colnames(AAOGrid_spec) <- AAOGrid_lbd


# Test sur les matrices 

dim(AAOGrid_spec)

dim(AAOGrid_CATAID)

length(AAOGrid_lbd)

AAOGrid_spec[1:5, 1:10]

range(AAOGrid_spec)

sum(is.na(AAOGrid_spec))


# Sauvegarde 

save(AAOGrid_CATAID, file="AAOGrid_CATAID.RData")

save(AAOGrid_lbd, file="AAOGrid_lbd.RData")

save(AAOGrid_spec, file="AAOGrid_spec.RData")



#===============================================================
# Visualisation des étapes du prétraitement d'un spectre AAOmega
#===============================================================

i <- 1   # choisir le spectre à afficher

# Spectre brut
lambda_brut <- lambdas_obs[[i]]
flux_brut <- flux_AAO[[i]]

# Spectre corrigé redshift
lambda_dered <- lambdas_corr_AAO[[i]]
flux_dered <- flux_AAO[[i]] 

# Spectre resamplé
flux_resamp <- flux_AAO_resampled_approx[[i]]

# Spectre normalisé
flux_norm <- flux_AAOnorm[[i]]


pdf("Etapes_pretraitement_spectre.pdf", width=14, height=10)

par(mfrow=c(2,2),
    mar=c(4,4,3,1)) 


#---------------------------------
# 1) Spectre brut
#---------------------------------

plot(lambda_brut,
     flux_brut,
     type="l",
     col="black",
     lwd=1,
     xlab=expression(lambda~(Angstrom)),
     ylab="Flux",
     main="(a) Spectre brut") 


#---------------------------------
# 2) Spectre corrigé redshift
#---------------------------------

plot(lambda_dered,
     flux_dered,
     type="l",
     col="red",
     lwd=1,
     xlab=expression(lambda~(Angstrom)),
     ylab="Flux",
     main="(b) Spectre corrigé du redshift")


#---------------------------------
# 3) Spectre resamplé
#---------------------------------

plot(lambda_AAOgrid,
     flux_resamp,
     type="l",
     col="blue",
     lwd=1,
     xlab=expression(lambda~(Angstrom)),
     ylab="Flux",
     main="(c) Spectre rééchantillonné")


#---------------------------------
# 4) Spectre normalisé
#---------------------------------

plot(lambda_AAOgrid,
     flux_norm,
     type="l",
     col="darkgreen",
     lwd=1,
     xlab=expression(lambda~(Angstrom)),
     ylab="Flux normalisé",
     main="(d) Spectre normalisé")


dev.off() 

