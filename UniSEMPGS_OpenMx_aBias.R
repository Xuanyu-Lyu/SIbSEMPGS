# the function script for fitting the OpenMx model with different levels of a
# adapted from Yongkong's function perform_SEM_model2_eq

fitDifferentA <- function(df, a, max.cores = 2){
    library(OpenMx)
    library(data.table)
    library(stringr)

  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")

  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'

  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")

  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  nv=1 #Number of variants
  dat_SEM <- df
  #dat[,c(NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel)]


  #Colnames
  colnames(dat_SEM) <- c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")
  nv=1

  #Scaling factor of the PRS - change depending on how PRS is scaled
  #k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS)=.5 at t0
  j <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_latent_hap_PRS_t0",name="j") #for var(hap_PRS_lat)=.5 at t0
  k <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=FALSE,values=.50,label="Var_hap_PRS_t0",name="k") #for var(hap_PRS_lat)=.5 at t0
  #k <- mxAlgebra(.5-2*g,name="k") #for var(hap_PRS) if the full PRS is scaled at equilibrium to have var(PRS) = 1


  # Paths for AFE model
  f <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.18,label="F",name="f",lbound=-.01)
  e <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.6,label="EnvPath",name="e",lbound=-.01)
  g <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=.01,label="hap_PRS_cov",name="g")
  delta <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=sqrt(.16),label="obs_coef",name="delta",lbound=-.01)
  a <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=FALSE,values=a,label="latent_coef",name="a",lbound=-.01)


  #Constraints
  #Matrices for latent variables & nonlinearly constrained estimates
  x1 <- mxMatrix(type="Full",nrow=nv,ncol=nv,free=TRUE,values=0.1,label="LatentF1",name="x1")
  w1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.15,label="covAandF",name="w1")
  mu1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.05,label="copath",name="mu1")
  i1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.02,label="covOL",name="i1")
  #sigma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=var(dat_SEM$Yo),label="sigma",name="sigma1")
  Omega1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.3,label="Omega",name="Omega1")
  Gamma1 <- mxMatrix(type="Lower",nrow=nv,ncol=nv,free=TRUE,values=0.5, label="Gamma",name="Gamma1")


  # mxAlgebra - nonlinear constraints
  x2 <- mxAlgebra(2*f^2*sigma1*(1 + sigma1*mu1),name="x2")
  w2 <- mxAlgebra(2*f*Omega1*(1 + sigma1*mu1),name="w2")
  v1 <- mxAlgebra((w1*a)/delta,name="v1")
  i2 <- mxAlgebra(Gamma1*mu1*Omega1,name="i2")
  g2 <- mxAlgebra(Omega1^2*mu1,name="g2")
  Omega2 <- mxAlgebra(2*a*i1 + delta*k + 2*delta*g+.5*w1,name="Omega2") #more general: between parent Y and their PRS
  Gamma2 <- mxAlgebra(2*a*h + a*j + 2*delta*i1 + .5*v1,name="Gamma2")
  #h <- mxAlgebra((g*a^2)/delta^2,name="h")
  h <- mxAlgebra(Gamma1^2*mu1,name="h")


  # Equating nonlinear constraints and parameters
  xCon <- mxConstraint(x1==x2,name="xCon")
  wCon <- mxConstraint(w1==w2,name="wCon")
  iCon <- mxConstraint(i1==i2,name="iCon") #not sure if needed - yes, seems to be needed
  gCon <- mxConstraint(g==g2,name="gCon")
  sigmaCon <- mxConstraint(sigma1==sigma2,name="sigmaCon")
  OmegaCon <- mxConstraint(Omega1==Omega2,name="OmegaCon")
  GammaCon <- mxConstraint(Gamma1==Gamma2,name="GammaCon")


  # mxAlgebra for implied parameters and relative covariances
  sigma1 <- mxAlgebra(2*a*Gamma1 + 2*delta*Omega1 + a*v1 + delta*w1 + x1 + e^2,name="sigma1")
  ThetaNTM05 <- mxAlgebra((Omega1-delta*k),name="ThetaNTM05") #this is HALF theta_NTM in the math
  spsPRS <-mxAlgebra(sigma1*mu1*Omega1,name="spsPRS") #THIS HAS BEEN FIXED; more general: between parent Y and spouse's PRS
  PO <- mxAlgebra((a*Gamma1 + delta*Omega1 + f*sigma1)*(1+mu1*sigma1),name="PO")
  spouse <-mxAlgebra(mu1*sigma1^2,name="spouse")


  #Covariance and Means Matrices
  CVmat<-    mxAlgebra(rbind(
    #       NTm       Tm          NTp     Tp          Ym         Yp        Yo
    cbind(  k+g        ,g         ,g          ,g        ,Omega1     ,spsPRS    ,ThetaNTM05   ),    #NTm
    cbind(  g          ,k+g       ,g          ,g        ,Omega1     ,spsPRS    ,Omega1        ),    #Tm
    cbind(  g          ,g         ,k+g        ,g        ,spsPRS     ,Omega1    ,ThetaNTM05   ),    #NTp
    cbind(  g          ,g         ,g          ,k+g      ,spsPRS     ,Omega1    ,Omega1     ),    #Tp
    cbind(  Omega1     ,Omega1    ,spsPRS     ,spsPRS   ,sigma1     ,spouse    ,PO         ),    #Ym
    cbind(  spsPRS     ,spsPRS    ,Omega1     ,Omega1   ,spouse     ,sigma1    ,PO         ),    #Yp
    cbind(  ThetaNTM05 ,Omega1    ,ThetaNTM05 ,Omega1   ,PO         ,PO        ,sigma1     ) ),  #Yo
    dimnames=list(colnames(dat_SEM),colnames(dat_SEM)),name="expCov")

  MNmat <- mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE,  values= rep(.01,nv), label=c("meanNTm","meanTm","meanNTp","meanTp","meanYm","meanYp","meanYo"),dimnames=list(NULL,c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")), name="expMean")


  #Objective object & parameters
  funML <- mxFitFunctionML()
  objdat <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("NTm","Tm","NTp","Tp","Ym","Yp","Yo"))

  params <- list(k ,j,  f , e , g ,delta,a,
                 x1 , w1 ,mu1 , i1, sigma1, Omega1 , Gamma1,
                 x2, w2, v1, i2, g2, Omega2 , Gamma2,
                 xCon, wCon, iCon, gCon, OmegaCon, GammaCon,
                 ThetaNTM05, h,
                 spsPRS, PO, spouse,
                 CVmat, MNmat,
                 funML,objdat)


  #Run the model
  dat_SEM2 <- mxData(observed=dat_SEM, type="raw")
  modelAFE <- mxModel("AFEmodel",params,dat_SEM2)
  AFE.Fit=try(mxRun(modelAFE,intervals=FALSE,silent=TRUE),silent=TRUE)
  if(class(AFE.Fit)!="try-error"){
    VAO=mxEval(delta^2*2*(k+2*g), AFE.Fit$AFEmodel, T)
    VAL=mxEval(a^2*2*(j+2*h), AFE.Fit$AFEmodel, T)
    VF=mxEval(x1, AFE.Fit$AFEmodel, T)
    VE=mxEval(e%*%e, AFE.Fit$AFEmodel, T)
    sigma=mxEval(sigma1,AFE.Fit$AFEmodel,T)
    west=mxEval(w1,AFE.Fit$AFEmodel,T)
    vest=mxEval(v1,AFE.Fit$AFEmodel,T)
    gest=mxEval(g2,AFE.Fit$AFEmodel,T)
    hest=mxEval(h,AFE.Fit$AFEmodel,T)
    h2=(VAO+VAL)/sigma
    VFratio=VF/sigma
    k_est=mxEval(k,AFE.Fit$AFEmodel,T)
    mu_est=mxEval(mu1,AFE.Fit$AFEmodel,T)
    f_est=mxEval(f,AFE.Fit$AFEmodel,T)
    deltaest=mxEval(delta,AFE.Fit$AFEmodel,T)
    aest=mxEval(a,AFE.Fit$AFEmodel,T)
    eest=mxEval(e,AFE.Fit$AFEmodel,T)
    ll=summary(AFE.Fit)$Minus2LogLikelihood
    results44=data.frame(method="model2_eq",VAO=VAO,VAL=VAL,VF=VF,VE=VE,sigma=sigma,west=west,vest=vest,gest=gest,hest=hest, h2=h2,k_est=k_est,mu=mu_est,f=f_est,
                         deltaest=deltaest,aest=aest,eest=eest,ll=ll)
  }else{
    results44=data.frame(method="model2_eq",VAO=NA,VAL=NA,VF=NA,VE=NA,sigma=NA,west=NA,vest=NA,gest=NA,hest=NA,h2=NA,k_est=NA,mu=NA,f=NA,
                         deltaest=NA,aest=NA,eest=NA,ll=NA)
  }

}
