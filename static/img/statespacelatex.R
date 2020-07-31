library(ctsem)
coupled <- ctModel(LAMBDA=matrix(c(1,0,0,1),2,2),DRIFT=c(-1,0,.5,-1),MANIFESTVAR=c(.1,0,0,.1),
  DIFFUSION=c(1,0,0,1),type='stanct',TDpredNames = 'GreatSpeech',latentNames = c('Motivation','Fitness'),
  TDPREDEFFECT=matrix(c(1,0)),manifestNames = c('ObsMotivation','ObsFitness'))
coupled$pars$indvarying<-FALSE
l=ctModelLatex(coupled,matrixnames = F,linearise = T)


indicators <- ctModel(LAMBDA=matrix(c(1,.97),2,1),DRIFT=c(-1),MANIFESTVAR=c(.1,0,0,.1),
  DIFFUSION=c(1),MANIFESTMEANS=c(0,1))
l=ctModelLatex(indicators,matrixnames = F,linearise = T)
cat(l)


common <- ctModel(LAMBDA=matrix(c(1,0,0,1),2,2),DRIFT=c(-1,0,0,-1),MANIFESTVAR=c(.1,0,0,.1),
  DIFFUSION=c(1,1,0,0))
l=ctModelLatex(common,matrixnames = F,linearise = T,compile = F)
cat(l)
