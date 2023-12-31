#-------------------------
#predict score for rt
#----------------------------
#' This function is used to predict the retention time of compounds and calculate scores
#' An excel document with positive controls is requried to build the calibration curves
#' @param Target 
#' @param Database 
#' @param RT.coeff 
#' @import rcdk
#' @return
#' @export
#'

ScoreRT<-function(Target,Database){
  
  #Using positive controls to build retention time predictions
  my_data <- read_excel("RTpred.xlsx")
  my_data$XlogP<-rep(0,nrow(my_data))
  for (i in 1:nrow(my_data)){
    mol<-my_data$SMILES[i]
    mol<-parse.smiles(mol)[[1]]
    convert.implicit.to.explicit(mol)
    my_data$XlogP[i]<-get.xlogp(mol)
  }
  
  #'regression
  results<-lm(RT~XlogP,data=my_data)
    S.single<-summary(results)
  RT.coeff<-S.single$coefficients
  SD.rt<-sd(S.single$residuals)
  
  #'predicting retention times of other compounds
  allscore0<-rep(0,nrow(Target))
  Target$rtscore<-rep(0,nrow(Target))
  for (i in 1:nrow(Target)){
    score<-Target$allscore[i]
    if (score==0){next}
    temp<-strsplit(as.character(score),';')
    temp<-as.numeric(unlist(temp))
    LogP<-NULL
    Tpsa<-NULL
    molsave<-unlist(strsplit(Target$SMILES[i],';'))
    for (j in 1:length(temp)){
      mol<-molsave[j]
      mol<-parse.smiles(mol)[[1]]
      convert.implicit.to.explicit(mol)
      LogP<-c(LogP,get.xlogp(mol))
      Tpsa<-c(Tpsa,get.tpsa(mol))
      predictRT<-RT.coeff[1]+LogP*RT.coeff[2]
      Diff.rt<-predictRT-Target$rt[i]
      temp.score<-log(exp(-0.25*Diff.rt^2/SD.rt^2))
    }
    Target$rtscore[i]<-paste(temp.score,collapse=';')
    allscore0[i]<-paste(temp+temp.score,collapse=';')
  }
  Target$allscore<-allscore0
  return(Target)
}