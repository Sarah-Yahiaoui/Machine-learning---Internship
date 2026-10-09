  library(FisherEM)
  Sys.info()["nodename"]
  sessionInfo()
  source("C:/Users/USER/Documents/mspsplit.R")
  
  #load("GAMA_SDSS_spec.RData")
  #load("GAMA_SDSS_lbd.RData")
  load("AAOmega_spectra_scaled.RData")
  load("AAOmega_lambda.RData")
   X <- X_scaled
   Y <- lambda
  length(lambda)
  length(Y)
  # /home/jmoultaka/Sarah/mspsplit.R, le chemin pour que les sorties vienne ici
  
  ## INPUT
  path = "./"  #/home/JIHANE/:wq Input/"
  model = "all" #c("AkjBk", "AkBk", "AjBk", "DkBk", "DBk", "ABk")
  K = 2
  method = "svd"
  maxit = 60
  nstart = 25
  nretry = 2 #5
  init = "random" #"kmeans"
  kernel = "linear"

  ## FUNCTIONS
  applyFEM <- function(spectra, lambdas, filename, method, maxit, nstart){
    ICL = -19999999 #300000
    isSolution = FALSE
    i=1
    while(i<=nretry){ #Boucle jusqu'à solution ou nretry essais
      #Application de FEM
      message(paste0("[",Sys.time(),"]\tApplying FEM: iteration ",i,"...\n"))
      try(assign(filename,fem(spectra,K=K,model=model,method=method,maxit=maxit,nstart=nstart,disp=T, crit="icl", init = init, kernel=kernel)))
      
      if(get(filename)$critValue==-Inf){ #Si aucune solution n'est trouvee (i.e. ICL==-Inf)
        i=i+1
        message("No solution found.\n")
        if(i<=nretry){message("Retrying...\n")} #Si le nb de tentatives n'a pas dépassé la limite
        else{message("Max iterations reached. Quiting without a solution...\n")} #Si le nb de tentatives a atteint la limite
      }
      else{ #Si une solution est trouvee (i.e. ICL!=-Inf)
	i = i + 1
        if (ICL < get(filename)$critValue) {
		      ICL = get(filename)$critValue
       		isSolution = TRUE
        	message("Solution found!\nK=",get(filename)$K,
                	",\tModel=",get(filename)$model, 
                	",\tICL=",get(filename)$critValue,"\nSaving...\n")
        
        	#try(save(list=filename,file=paste0("/home/jmoultaka/Results_FEM/",filename,".RData")))
        	#pdf(file=paste0("/home/jmoultaka/Results_FEM/",filename,".pdf"))
        	try(save(list=filename,file=paste0("./",filename,".RData")))
                pdf(file=paste0("./",filename,".pdf"))
		try(mspsplit(get(filename),spec=spectra,Xlambdas=lambdas,yrange=c(0.2,2),legcex=0.5,gpord=as.numeric(names(sort(table(get(filename)$cls),decreasing = T)))))
        	try(plot(get(filename)))
        	dev.off()
	}
      }
      rm(list=filename)
    }
    message("---------------------------------------------------\n") #Fin du bin
  }
  is_same_length <- function(files_lbd, files_spec){
    return(length(files_lbd)==length(files_spec))
  }
  is_same_names <- function(file_lbd, file_spec){
    name_lbd = gsub("_lbd.RData","", file_lbd)
    name_spec = gsub("_spec.RData","", file_spec)
    return(name_lbd==name_spec)
  }
  get_name <- function(file){
    res = gsub("_spec.RData","", file)
    res = gsub("_lbd.RData","", res)
    return(res)
  }
  
  ## MAIN
  files_lbd = list.files(path, pattern="lbd.RData$") ; files_spec = list.files(path, pattern="spec.RData$")
  #files_lbd = rev(files_lbd) ; files_spec = rev(files_spec) # inversement de la liste des fichiers pour commencer aux z >>
  ## MAIN
  message("\n---------------------------------------------------\n")
  message(paste0("Model=",model))
  message(paste0("K=",K[1],":",K[length(K)]))
  message(paste0("method=",method))
  message(paste0("maxit=",maxit))
  message(paste0("nstart=",nstart))
  message(paste0("nretry=",nretry))
  message(paste0("init=",init))
  

 

  if (!is_same_length(files_lbd, files_spec)){message("/!\\ \tErreur à l'initialisation des données : fichiers manquants.")
  }else{
    for (i in 1:length(files_lbd)){
      file_lbd = files_lbd[i] ; file_spec = files_spec[i]
      if (!is_same_names(file_lbd, file_spec)){message(" /!\\ \t Erreur à l initialisation des données : noms de fichiers de lbd et spec différents.\n")
      }else{
        message("[",Sys.time(),"]\tChargement de ",get_name(file_lbd))
        spec = X  ; rm(list=gsub(".RData","",file_spec)) ; lbd = Y ; rm(list=gsub(".RData","",file_lbd))
        message("[",Sys.time(),"]\tN spectres : ", nrow(spec),"\tN longueurs d'ondes : ", ncol(spec))
        applyFEM(spectra=spec, lambdas=lbd, filename=get_name(file=file_lbd),method=method, maxit=maxit, nstart=nstart)
        rm(spec) ; rm(lbd)
      }
    }
  } 

  

  