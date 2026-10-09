# le code pour avoir le fichier pdf pour les clusters des spectres 
# dispersion = quantiles or std

mspsplit <- function(FEMout,spec=spectra,Xlambdas=lambdas,gpord=c(1:FEMout$K),arrang=NULL,disp="quantiles",yrange=NULL,legcex=0.8,legwidth=NULL) {

  # La partie ou je recupere le spectre moyen
  
    if(length(FEMout$my)==0) specmeans <- FEMout$mean else specmeans <- FEMout$my #  if version==151 else version=14  
  #    K <- FEMout$K
  ynull <- F
  if (is.null(yrange)) ynull=T
  #if (is.null(arrang)) arrang <- c(ceiling(sqrt(length(gpord))),floor(sqrt(length(gpord))))  
  
  if (is.null(arrang)) {
    nK <- min(length(gpord),length(table(FEMout$cls)))
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if(is.wholenumber(sqrt(nK))) arrang <- c(sqrt(nK),sqrt(nK)) else {
      a <- ceiling(sqrt(nK))
      arrang <- c(a,ceiling(nK/a))
    }
  
      # {
    # cof <- function(a,K) a*a+a-K
    # tmp <- function(a,K) min(a[which(min(cof(a,K)) & cof(a,K)>=0)])
    # 
    # arrang <- c(tmp(1:nK,nK)+1,tmp(1:nK,nK))      }
  }
  
  library(FisherEM)
  library(scales) # pour la fonction alpha() de la transparence
  #  op <- par(mfrow=arrang,mar=c(3,3,1,1),mgp=c(1.5,0.5,0),oma=c(3,0,2,0)) # ou mar=c(3,3,0.5,0.5),mgp=c(1.5,0,0)
  op <- par(mfrow=arrang,mar=c(0,0,0,0),mgp=c(1.5,0.5,0),oma=c(5,4,1,1),cex=1)
  j <- 0
  for (icls in gpord) {
    i <- which(as.numeric(names(table(FEMout$cls)))==icls) ; if (length(i)==0) {next}
    cof <- cbind(Xlambdas,specmeans[icls,])
    colnames(cof) = c("lambda","flux") 
    #      plot(cof,type="l",ylim=yrange)    
    if (table(FEMout$cls)[i] != 1) { # pour ??liminer les groupes avec un seul objet
      cof1 <- as.matrix(spec)[which(FEMout$cls==icls),]
      
      if (disp == "quantiles") {
        stdmin <- cbind(Xlambdas,apply(cof1,2,function(x) quantile(x,probs=0.1)))
        stdmax <- cbind(Xlambdas,apply(cof1,2,function(x) quantile(x,probs=0.9))) 
      }
      if (disp == "std") {
        cof1std <- cbind(Xlambdas,apply(cof1,2,sd)) 
        stdmin <- cof ; stdmin[,2] <- cof[,2]-cof1std[,2]
        stdmax <- cof ; stdmax[,2] <- cof[,2]+cof1std[,2]
      }
      colnames(stdmin) = c("lambda","flux")
      colnames(stdmax) = c("lambda","flux")
      if (ynull) yrange=c(min(stdmin[,2]),max(stdmax[,2]))
      #      plot(cof,type="l",ylim=yrange)    
      plot(cof,type="l",ylim=yrange,xlab="",ylab="",axes=F)    
      lines(stdmin,col=alpha("gray",0.6)) 
      lines(stdmax,col=alpha("gray",0.6)) 
      #title(main=paste0("Group ",i))
      if (is.null(legwidth)) legwidth=min(1,2000/strwidth("mean G N  N=100000 "))
      if (disp == "quantiles") {legend("topleft",legend=c(paste("mean G",icls,"  N=",table(FEMout$cls)[i]),"10% and 90% quantiles"),col=c("black",alpha("gray",0.6)),lty=1,bty="n",cex=legwidth) }
      if (disp == "std") {legend("topleft",legend=c(paste("mean G",icls,"  N=",table(FEMout$cls)[i]),"Standard deviation"),col=c("black",alpha("gray",0.6)),lty=1,bty="n",cex=legwidth) }
      polygon(c(Xlambdas,rev(Xlambdas)),c(stdmax[,2],rev(stdmin[,2])),col=alpha("gray",0.3),border=NA,ylim=c(0.3,2))
    }
    if (table(FEMout$cls)[i] == 1) {       
      #     plot(cof,type="l",ylim=yrange)   
      plot(cof,type="l",ylim=yrange,xlab="",ylab="",axes=F)   
      if (is.null(legwidth)) legwidth=min(1,2000/strwidth("mean G N  N=100000 "))
      legend("topleft",legend=paste("mean G",icls,"  N=",table(FEMout$cls)[i]),col="black",lty=1,bty="n",cex=legwidth)}
    box(which = "plot", bty = "o") 
    if (j +1 > length(gpord)-arrang[2])  axis(side = 1,cex.axis=legcex) 
    if (j%%arrang[2] == 0) axis(side = 2,cex.axis=legcex)
    j <- j+1
  }
  par(op)
}

