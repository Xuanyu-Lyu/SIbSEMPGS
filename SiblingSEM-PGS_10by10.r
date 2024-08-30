# This is a OpenMx script to fit a SEM-PGS model with parental and sibling data
# This is a trial version of 10 by 10 matrix due to the linear dependency issue
# Author: Xuanyu Lyu
# Date: 05/20/2024


# load the packages
library(OpenMx)
library(data.table)
library(stringr)
#source('http://openmx.ssri.psu.edu/getOpenMx.R')

# Specify Options:
    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")

# Load the simulated data for this demonstration:
    Example_Data  <- fread("Example_Data.txt", header = T)
    Example_Data <- Example_Data[,-c(10,12)]
#    str(Example_Data)
 cov(Example_Data, use="pairwise.complete.obs")
# Create variables-- For some, we will also input their algebraic expectations, which we can obtain from path tracing
    # Variance Components:

    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0, label="VY1", name="VY") # Phenotypic variance
	VF    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0, label="VF1", name="VF") # Variance due to VT
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0, label="VE1", name="VE") # Residual variance

    VY_Algebra <- mxAlgebra(2 * Omega * delta + delta * w + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY * (1 + VY * mu),            name="VF_Algebra")

    # Genetic effects:
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="delta1",name="delta") # Effect of PGS on phen
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k1",    name="k")     # PGS variance (if no AM)
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.3, label="Omega1",name="Omega") # Within-person PGS-Phen covariance

    Omega_Algebra <- mxAlgebra(2 * delta * g + delta * k + .5 * w, name="Omega_Algebra") # E.g., cov(Yp, NTp)

    # Assortative mating effects:
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.4, label="mu1", name="mu") # AM co-path coefficient
    g     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.05, label="g1",  name="g")  # Increase in cross-mate PGS (co)variances from AM

    # Vertical transmission effects (note that VF above is also due to VT):
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=0, label="f1",  name="f") # Vertical Transmission effect
    w     <- mxAlgebra(2 * f * Omega * (1 + VY*mu), name="w")                                # Genetic nurture

# Set the variables created above to be equal to their algebraic constraints:
    VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
    VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
	Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')
    Phen_Homog <- mxConstraint(g == Omega^2 * mu,  name="Phen_Homog")

# Provide the algebraic expectations for the covariances between relatives:
    # Covariances between the offspring 1 phenotype and offspring 2 phenotype:
    Yo1_Yo2 <- mxAlgebra(2 * delta * w + 4 * delta^2 * g + delta^2 * k + VF, name="Yo1_Yo2")
    # Covariances between the offspring 1 phenotype and parental PGS's:
    Yo1_NTp1 <- mxAlgebra(2 * delta * g + .5 * w, name="Yo1_NTp1")
    Yo1_NTm1 <- mxAlgebra(2 * delta * g + .5 * w, name="Yo1_NTm1")
    Yo1_Tp1  <- mxAlgebra(Yo1_NTp1 + delta * k,     name="Yo1_Tp1")
    Yo1_Tm1  <- mxAlgebra(Yo1_NTm1 + delta * k,     name="Yo1_Tm1")
    #Yo1_NTp2 <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w, name="Yo1_NTp2")
    #Yo1_NTm2 <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w, name="Yo1_NTm2")
    Yo1_Tp2  <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w, name="Yo1_Tp2")
    Yo1_Tm2  <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w,name="Yo1_Tm2")

    # Covariances between the offspring 1 phenotype and parental phenotypes:
    Yo1_Ym <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo1_Ym")
    Yo1_Yp <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo1_Yp")

    # Covariances between the offspring 2 phenotype and parental PGS's:
    #Yo2_NTp2 <- mxAlgebra(2 * delta * g + .5 * w, name="Yo2_NTp2")
    #Yo2_NTm2 <- mxAlgebra(2 * delta * g + .5 * w, name="Yo2_NTm2")
    Yo2_Tp2  <- mxAlgebra(2 * delta * g + .5 * w + delta * k, name="Yo2_Tp2")
    Yo2_Tm2  <- mxAlgebra(2 * delta * g + .5 * w + delta * k, name="Yo2_Tm2")
    Yo2_NTp1 <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w, name="Yo2_NTp1")
    Yo2_NTm1 <- mxAlgebra(2 * delta * g + .5*delta*k + .5 * w, name="Yo2_NTm1")
    Yo2_Tp1  <- mxAlgebra(Yo2_NTp1,     name="Yo2_Tp1")
    Yo2_Tm1  <- mxAlgebra(Yo2_NTm1,     name="Yo2_Tm1")

    # Covariances between the offspring 2 phenotype and parental phenotypes:
    Yo2_Ym <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo2_Ym")
    Yo2_Yp <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo2_Yp")

    # Between-spouse covariances:
    Yp_PGSm <- mxAlgebra(VY * mu * Omega, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(VY * mu * Omega, name="Ym_PGSp")
    Ym_Yp   <- mxAlgebra(VY * mu * VY,    name="Ym_Yp")



# Expected covariances between each variable:
    CovMatrix <-    mxAlgebra(rbind(
    #       Yo1        Yo2       Yp        Ym        Tp1       NTp1      Tm1         NTm1       Tp2           Tm2            
    cbind(  VY      ,  Yo1_Yo2  ,Yo1_Yp   ,Yo1_Ym   ,Yo1_Tp1  ,Yo1_NTp1 ,Yo1_Tm1     ,Yo1_NTm1 ,Yo1_Tp2    ,Yo1_Tm2       ),    #Yo1
    cbind(  Yo1_Yo2 ,  VY       ,Yo2_Yp   ,Yo2_Ym   ,Yo2_Tp1  ,Yo2_NTp1 ,Yo2_Tm1     ,Yo2_NTm1 ,Yo2_Tp2    ,Yo2_Tm2       ),    #Yo2
    cbind(  Yo1_Yp  ,  Yo2_Yp   ,VY       ,Ym_Yp    ,Omega    ,Omega    ,Yp_PGSm     ,Yp_PGSm  ,Omega      ,Yp_PGSm       ),    #Yp
    cbind(  Yo1_Ym  ,  Yo2_Ym   ,Ym_Yp    ,VY       ,Ym_PGSp  ,Ym_PGSp  ,Omega       ,Omega    ,Ym_PGSp    ,Omega         ),    #Ym
    cbind(  Yo1_Tp1 ,  Yo2_Tp1  ,Omega    ,Ym_PGSp  ,k+g      ,g        ,g           ,g        ,.5*k+g     ,g             ),    #Tp1
    cbind(  Yo1_NTp1,  Yo2_NTp1 ,Omega    ,Ym_PGSp  ,g        ,k+g      ,g           ,g        ,.5*k+g     ,     g             ),    #NTp1
    cbind(  Yo1_Tm1 ,  Yo2_Tm1  ,Yp_PGSm  ,Omega    ,g        ,g        ,k+g         ,g        ,g          ,.5*k+g             ),    #Tm1
    cbind(  Yo1_NTm1,  Yo2_NTm1 ,Yp_PGSm  ,Omega    ,g        ,g        ,g           ,k+g      ,g          ,.5*k+g             ),    #NTm1
    cbind(  Yo1_Tp2 ,  Yo2_Tp2  ,Omega    ,Ym_PGSp  ,.5*k+g  ,.5*k+g  ,  g     ,      g        ,k+g        ,g                  ),    #Tp2
    #cbind(  Yo1_NTp2,  Yo2_NTp2 ,Omega    ,Ym_PGSp  ,.5*k+g ,.5*k+g  ,  g     ,      g        ,g          ,g                  ),    #NTp2
    cbind(  Yo1_Tm2,   Yo2_Tm2  ,Yp_PGSm  ,Omega    ,     g  , g        ,.5*k+g     ,.5*k+g    ,g          ,k+g                )),    #Tm2
    #cbind(  Yo1_NTm2,  Yo2_NTm2 ,Yp_PGSm  ,Omega    ,     g , g        ,.5*k+g     ,.5*k+g    ,g          ,g                 )),   #NTm2
    dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

# Expected means for each variable:
    Means <-   mxMatrix(type="Full", nrow=1, ncol=10, free=TRUE, values= .1, label=c("meanYo1","meanYo2","meanYp","meanYm","meanTp1","meanNTp1","meanTm1","meanNTm1","meanTp2","meanTm2"), 
                        dimnames=list(NULL, c("Yo1","Yo2","Yp","Ym","Tp1","NTp1","Tm1","NTm1","Tp2","Tm2")), name="expMean")

# Connect the variable means and covariances with one another:
    ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("Yo1","Yo2","Yp","Ym","Tp1","NTp1","Tm1","NTm1","Tp2","Tm2"))

# Convert data into a usable format for OpenMx:
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

# Create fit function:
    FitFunctionML <- mxFitFunctionML()

# Specify what parameters we're going to be including in our model:
Params <- list( VY, VY_Algebra, VY_Constraint, VF, VF_Algebra, VF_Constraint, VE,  # Variance Components
                Omega, Omega_Algebra, Omega_Constraint, delta, k,                    # Genetic Effects
                mu, g,                                                              # AM Effects
                f, w,                                                                # VT Effects
                Yo1_Yp, Yo1_Ym, Yo1_Tp1, Yo1_Tm1, Yo1_NTp1, Yo1_NTm1, Yo1_Tp2, Yo1_Tm2, # Relative Covariances for Offspring 1
                Yo2_Yp, Yo2_Ym, Yo2_Tp1, Yo2_Tm1, Yo2_NTp1, Yo2_NTm1, Yo2_Tp2, Yo2_Tm2,  # Relative Covariances for Offspring 2
                Yo1_Yo2, Yp_PGSm, Ym_PGSp, Ym_Yp, # Relative Covariances
                FitFunctionML, Means, ModelExpectations, CovMatrix)                  # Model Parameters

# Create the model:
    options(warning.length = 8000)
    Model1 <- mxModel("SiblingSEM_PGS", Params, Example_Data_Mx, Phen_Homog)

    fitModel1 <- mxRun(Model1,   intervals=F, silent=F)


cov(Example_Data, use="pairwise.complete.obs")
