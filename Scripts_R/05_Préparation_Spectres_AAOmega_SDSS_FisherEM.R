
#==================================================================================
# les AAOmega et sdss sur la meme grille 
#====================================================================================

load("AAOsdssGrid_CATAID.RData") 

load("AAOsdssGrid_lbd.RData")

load("AAOsdssGrid_spec.RData")


# Construction de la matrice SDSS 

#Get file names from directory

directory <- "C:/Users/USER/Desktop/Sarah/GAMA_SDSS_Normalized"

files <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
length(files) # 2550

#Split to save names ; name for data frame will be first element

names <- strsplit(files, "\\.txt$")


#Matrix set-up

GAMA_SDSS_spec <- matrix(nrow = length(names), ncol = 3681) #NB of spectra x NB of wavelengths
GAMA_SDSS_CATAID <- matrix(nrow = 2550, ncol = 1) #CATAID
dim(GAMA_SDSS_CATAID)# 2550    1
dim(GAMA_SDSS_spec) # 2550 3681

for (i in 1:length(files)) { #For each file in the list
  
  fileName <- files[[i]] #Save filename of element i
  
  tempData <- read.table(file = fileName, check.names = FALSE, header = TRUE) #Read txt file
  
  GAMA_SDSS_spec[i,] <- tempData[,1]
  
  GAMA_SDSS_CATAID[i] <- colnames(tempData)
  
  colnames(GAMA_SDSS_spec) <- rownames(tempData)
}

rownames(GAMA_SDSS_spec) <- GAMA_SDSS_CATAID[,1]

GAMA_SDSS_lbd <- as.numeric(rownames(tempData))

dim(GAMA_SDSS_CATAID)
dim(GAMA_SDSS_spec)
length(GAMA_SDSS_lbd)

length(intersect(AAOsdssGrid_CATAID, GAMA_SDSS_CATAID)) #276 

all.equal(GAMA_SDSS_lbd, AAOsdssGrid_lbd)


AAO_id  <- as.character(AAOsdssGrid_CATAID)
SDSS_id <- as.character(GAMA_SDSS_CATAID)

# Galaxies communes
commun <- intersect(AAO_id, SDSS_id)

length(commun) # 276

# Garder uniquement les AAOmega non présents dans SDSS

idx_AAO <- !(AAO_id %in% commun)

AAO_spec_unique <- AAOsdssGrid_spec[idx_AAO, ]
AAO_id_unique   <- AAO_id[idx_AAO]

# Vérification
length(intersect(AAO_id_unique, SDSS_id))


# la base finale 

AAO_SDSS_spec <- rbind(
  GAMA_SDSS_spec,
  AAO_spec_unique
)

AAO_SDSS_CATAID <- c(
  SDSS_id,
  AAO_id_unique
)

AAO_SDSS_lbd <- GAMA_SDSS_lbd

dim(AAO_SDSS_spec)
length(AAO_SDSS_CATAID)
length(AAO_SDSS_lbd)

nrow(AAO_SDSS_spec) == length(AAO_SDSS_CATAID)

AAO_SDSS_CATAID <- as.matrix(AAO_SDSS_CATAID)
dim(AAO_SDSS_CATAID)

range(AAO_SDSS_spec)
summary(as.vector(AAO_SDSS_spec)) 

# Sauvegarde des matrices

save(AAO_SDSS_CATAID, file="AAO_SDSS_CATAID.RData")

save(AAO_SDSS_spec, file="AAO_SDSS_spec.RData")

save(AAO_SDSS_lbd, file="AAO_SDSS_lbd.RData")




