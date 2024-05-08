
# Estimating Vertical Transmission and Genetic Nurture using Polygenic Scores #
# Jared V. Balbona, Yongkang Kim, and Matthew C. Keller #

# If you have any questions, please feel free to email Yongkang and/ or Jared at jared.balbona@colorado.edu and yongkangkim87@gmail.com.
# And for further information, please also check out our recent paper (doi.org/10.1007/s10519-020-10032-w) and the OpenMx website (https://openmx.ssri.psu.edu/).

# Thank you for following our demonstration!

#########
## KEY ##
#########
# Our dataset contains:
    # Yo:  The offspring's phenotypic value
    # Yp:  The father's phenotypic value
    # Ym:  The mother's phenotypic value
    # Tp:  A PGS made from the TRANSMITTED portion of the PATERNAL genotype
    # Tm:  A PGS made from the TRANSMITTED portion of the MATERNAL genotype
    # NTp: A PGS made from the NON-TRANSMITTED portion of the PATERNAL genotype
    # NTm: A PGS made from the NON-TRANSMITTED portion of the MATERNAL genotype

# Our model estimates:
    # VY: Phenotypic variance
    # VA: Variance due to additive
    # VF: Variance due to parental effects/ vertical transmissino
    # VE: Variance due to the offspring's unique environment
    # w:  Genetic Nurture (i.e., the passive G-E covariance due to vertical transmission)
    # f:  The direct impact of the parental phenotype on the offspring phenotype (i.e., vertical transmission)
    # g:  The increase in haplotypic (co)variances due to assortative mating
    # mu: The assortative mating co-path coefficient, equal to (Spousal Phen Cov) / (Spousal Phen Variance)^2
    # k:  The haplotypic PGS variance under no AM (a constant that is not estimated)
    # delta: The direct effect of the PGS on the phenotype
    # Omega: The covariance between a parent's phenotype and either of their haplotypic PGS's (e.g., cov(Ym, NTm))
    # h2: The proportion of phenotypic variance due to additive genetic factors (i.e., narrow-sense h2)
    # ll: The log-likelihood for our model. It essentially measures how much unexplained variation there is in our model, with higher values indicating less accuracy.

############################################################################################
# STEP O: Load data, specify options, and look at the observed variance/ covariance matrix #
############################################################################################
library(OpenMx)
library(data.table)
library(stringr)

# omxGetNPSOL() # Installs the NPSOL optimizer, which is required to estimate non-linear algebraic constraints (i.e., constraints where two variables are functions of each other). You won't need to run this command today for our workshop practical, but you may need to run it if you try to use this scripts outside of this workshop.

# Specify Options:
    mxOption(NULL,"Calculate Hessian","Yes")
    mxOption(NULL,"Standard Errors","Yes")
    mxOption(NULL,"Default optimizer","NPSOL")

# Load the simulated data for this demonstration:
    Example_Data  <- fread("Example_Data.txt", header = T)
    str(Example_Data)

# Check the column names-- For this script to work as intended, our data must contain the headers Yo, Yp, Ym, Tp, NTp, Tm, and NTm, with these exact spellings/ capitalizations and in this exact order
    colnames(Example_Data) == c("Yo", "Yp", "Ym", "Tp", "NTp", "Tm", "NTm")

# Look at the full observed variance/ covariance matrix:
    round(cov(Example_Data),3)

##################################################
# STEP 1: Create Variables and their Constraints #
##################################################

# Create variables-- For some, we will also input their algebraic expectations, which we can obtain from path tracing
    # Variance Components:

    VY    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="VY1", name="VY") # Phenotypic variance
	VF    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="VF1", name="VF") # Variance due to VT
    VE    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="VE1", name="VE") # Residual variance

    VY_Algebra <- mxAlgebra(2 * Omega * delta + delta * w + VF + VE, name="VY_Algebra")
    VF_Algebra <- mxAlgebra(2 * f^2 * VY * (1 + VY * mu),            name="VF_Algebra")

    # Genetic effects:
    delta <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="delta1",name="delta") # Effect of PGS on phen
    k     <- mxMatrix(type="Full", nrow=1, ncol=1, free=F, values=.5, label="k1",    name="k")     # PGS variance (if no AM)
    Omega <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="Omega1",name="Omega") # Within-person PGS-Phen covariance

    Omega_Algebra <- mxAlgebra(2 * delta * g + delta * k + .5 * w, name="Omega_Algebra") # E.g., cov(Yp, NTp)

    # Assortative mating effects:
    mu    <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="mu1", name="mu") # AM co-path coefficient
    g     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="g1",  name="g")  # Increase in PGS (co)variances from AM

    # Vertical transmission effects (note that VF above is also due to VT):
    f     <- mxMatrix(type="Full", nrow=1, ncol=1, free=T, values=.2, label="f1",  name="f") # Vertical Transmission effect
    w     <- mxAlgebra(2 * f * Omega * (1 + VY*mu), name="w")                                # Genetic nurture

# Set the variables created above to be equal to their algebraic constraints:
    VY_Constraint    <- mxConstraint(VY == VY_Algebra,       name='VY_Constraint')
    VF_Constraint    <- mxConstraint(VF == VF_Algebra,       name='VF_Constraint')
	Omega_Constraint <- mxConstraint(Omega == Omega_Algebra, name='Omega_Constraint')

# Provide the algebraic expectations for the covariances between relatives:
    # Covariances between the offspring phenotype and parental PGS's:
    Yo_NTp <- mxAlgebra(2 * delta * g + .5 * w, name="Yo_NTp")
    Yo_NTm <- mxAlgebra(2 * delta * g + .5 * w, name="Yo_NTm")
    Yo_Tp  <- mxAlgebra(Yo_NTp + delta * k,     name="Yo_Tp")
    Yo_Tm  <- mxAlgebra(Yo_NTm + delta * k,     name="Yo_Tm")

    # Covariances between the offspring phenotype and parental phenotypes:
    Yo_Ym <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo_Ym")
    Yo_Yp <- mxAlgebra(delta * Omega + delta * Omega * mu * VY + f * VY + f * VY * mu * VY, name="Yo_Yp")

    # Between-spouse covariances:
    Yp_PGSm <- mxAlgebra(VY * mu * Omega, name="Yp_PGSm")
    Ym_PGSp <- mxAlgebra(VY * mu * Omega, name="Ym_PGSp")
    Ym_Yp   <- mxAlgebra(VY * mu * VY,    name="Ym_Yp")

# Constraints for testing different model assumptions:
	No_VT      <- mxConstraint(f == 0,  name="No_VT") # Inclusion of this constraint will result in a model that assumes no VT
    Phen_Homog <- mxConstraint(g == Omega^2 * mu,  name="Phen_Homog") # Inclusion of this assumes AM is based on phenotypic similarity

#######################################################################################
# STEP 2: Relate our algebraic expectations to the observed data and create the model #
#######################################################################################

# Expected covariances between each variable:
    CovMatrix <-    mxAlgebra(rbind(
    #       Yo       Yp       Ym       Tp        NTp      Tm         NTm
    cbind(  VY      ,Yo_Yp    ,Yo_Ym    ,Yo_Tp	  ,Yo_NTp   ,Yo_Tm     ,Yo_NTm  ),    #Yo
    cbind(  Yo_Yp   ,VY       ,Ym_Yp    ,Omega    ,Omega    ,Yp_PGSm   ,Yp_PGSm ),    #Yp
    cbind(  Yo_Ym   ,Ym_Yp    ,VY       ,Ym_PGSp  ,Ym_PGSp  ,Omega     ,Omega   ),    #Ym
    cbind(  Yo_Tp   ,Omega    ,Ym_PGSp  ,k+g      ,g        ,g         ,g       ),    #Tp
    cbind(  Yo_NTp  ,Omega    ,Ym_PGSp  ,g        ,k+g      ,g         ,g       ),    #NTp
    cbind(  Yo_Tm   ,Yp_PGSm  ,Omega    ,g        ,g        ,k+g       ,g       ),    #Tm
    cbind(  Yo_NTm  ,Yp_PGSm  ,Omega    ,g        ,g        ,g         ,k+g     )),   #NTm
    dimnames=list(colnames(Example_Data),colnames(Example_Data)),name="expCov")

# Expected means for each variable:
    Means <-   mxMatrix(type="Full", nrow=1, ncol=7, free=TRUE, values= .1, label=c("meanYo","meanYp","meanYm","meanTp","meanNTp","meanTm","meanNTm"), dimnames=list(NULL, c("Yo","Yp","Ym","Tp","NTp","Tm","NTm")), name="expMean")

# Connect the variable means and covariances with one another:
    ModelExpectations <- mxExpectationNormal(covariance="expCov",means="expMean", dimnames=c("Yo","Yp","Ym","Tp","NTp","Tm", "NTm"))

# Convert data into a usable format for OpenMx:
    Example_Data_Mx <- mxData(observed=Example_Data, type="raw" )

# Create fit function:
    FitFunctionML <- mxFitFunctionML()

#####################################
# STEP 3: Run the model iteratively #
#####################################

Start.Val.List <- VT.AnyAM.Result <- VT.OnlyPH.Result <- NoVT.AnyAM.Result <- NoVT.OnlyPH.Result <- data.frame()
i = 1 # Used for indexing below

for(f.Start.Val in c(0, .4)){
  for(VE.Start.Val in c(.2, .5, .8)){
    for(delta.Start.Val in c(0, .5)){

      # Apply the starting values for this iteration:
      Start.Val.List = rbind(Start.Val.List, data.frame(f.Start.Val=f.Start.Val, VE.Start.Val=VE.Start.Val, delta.Start.Val=delta.Start.Val))

      VE    <- mxMatrix(type="Full", nrow=1	,ncol=1, free=TRUE, values=VE.Start.Val,    label="VE1",   name="VE")
      delta <- mxMatrix(type="Full", nrow=1	,ncol=1, free=TRUE, values=delta.Start.Val, label="delta1",name="delta")
      f     <- mxMatrix(type="Full", nrow=1	,ncol=1, free=TRUE, values=f.Start.Val,     label="f1",    name="f")

      # Specify what parameters we're going to be including in our model:
      Params <- list( VY, VY_Algebra, VY_Constraint, VF, VF_Algebra, VF_Constraint, VE,  # Variance Components
                    Omega, Omega_Algebra, Omega_Constraint, delta, k,                    # Genetic Effects
                    mu, g,                                                               # AM Effects
                    f, w,                                                                # VT Effects
                    Yo_Yp, Yo_Ym, Yo_Tp, Yo_Tm, Yo_NTp, Yo_NTm, Yp_PGSm, Ym_PGSp, Ym_Yp, # Relative Covariances
                    FitFunctionML, Means, ModelExpectations, CovMatrix)                  # Model Parameters

      # Construct each model based on their different assumptions; to keep track, name them based on the iteration number
	  assign(paste0("VT.AnyAM.It",i),   mxModel("VT_AnyAM",   Params, Example_Data_Mx                   )) # VT and g both freely estimable
	  assign(paste0("VT.OnlyPH.It",i),  mxModel("VT_OnlyPH",  Params, Example_Data_Mx, Phen_Homog       )) # Assumes g = Omega^2 mu
	  assign(paste0("NoVT.AnyAM.It",i), mxModel("NoVT_AnyAM", Params, Example_Data_Mx, No_VT            )) # Assumes No VT
	  assign(paste0("NoVT.OnlyPH.It",i),mxModel("NoVT_OnlyPH",Params, Example_Data_Mx, No_VT, Phen_Homog)) # Assumes No VT and g = Omega^2 mu

	  # Fit each of those different models to our data-- Again, we'll give the output a name that includes the iteration number.
	  assign(paste0("VT.AnyAM.It",i,".Fit"),    try(mxRun(eval(parse(text=paste0("VT.AnyAM.It",i))),    intervals=F, silent=T), silent=T))
	  assign(paste0("VT.OnlyPH.It",i,".Fit"),   try(mxRun(eval(parse(text=paste0("VT.OnlyPH.It",i))),   intervals=F, silent=T), silent=T))
	  assign(paste0("NoVT.AnyAM.It",i,".Fit"),  try(mxRun(eval(parse(text=paste0("NoVT.AnyAM.It",i))),  intervals=F, silent=T), silent=T))
	  assign(paste0("NoVT.OnlyPH.It",i,".Fit"), try(mxRun(eval(parse(text=paste0("NoVT.OnlyPH.It",i))), intervals=F, silent=T), silent=T))

	  # If the model ran without errors in this iteration, then record its estimated log likelihood; if not, then assign NA
	  VT.AnyAM.ll <- ifelse(class(eval(parse(text=paste0("VT.AnyAM.It",i,".Fit")))) != "try-error", summary(eval(parse(text=paste0("VT.AnyAM.It",    i, ".Fit"))))$Minus2LogLikelihood, NA)

	  VT.OnlyPH.ll   <- ifelse(class(eval(parse(text=paste0("VT.OnlyPH.It",  i,".Fit")))) != "try-error", summary(eval(parse(text=paste0("VT.OnlyPH.It",   i, ".Fit"))))$Minus2LogLikelihood, NA)

      NoVT.AnyAM.ll  <- ifelse(class(eval(parse(text=paste0("NoVT.AnyAM.It", i,".Fit")))) != "try-error", summary(eval(parse(text=paste0("NoVT.AnyAM.It",  i, ".Fit"))))$Minus2LogLikelihood, NA)

      NoVT.OnlyPH.ll    <- ifelse(class(eval(parse(text=paste0("NoVT.OnlyPH.It",   i,".Fit")))) != "try-error", summary(eval(parse(text=paste0("NoVT.OnlyPH.It", i, ".Fit"))))$Minus2LogLikelihood, NA)

      # Record the iteration number, log likelihood, and start values for the current iteration:
      VT.AnyAM.Result1    <- data.frame(rbind(c(Iteration=i, ll=VT.AnyAM.ll,    f.Start.Val=f.Start.Val, VE.Start.Val=VE.Start.Val, delta.Start.Val=delta.Start.Val)))
      VT.OnlyPH.Result1   <- data.frame(rbind(c(Iteration=i, ll=VT.OnlyPH.ll,   f.Start.Val=f.Start.Val, VE.Start.Val=VE.Start.Val, delta.Start.Val=delta.Start.Val)))
      NoVT.AnyAM.Result1  <- data.frame(rbind(c(Iteration=i, ll=NoVT.AnyAM.ll,  f.Start.Val=f.Start.Val, VE.Start.Val=VE.Start.Val, delta.Start.Val=delta.Start.Val)))
      NoVT.OnlyPH.Result1 <- data.frame(rbind(c(Iteration=i, ll=NoVT.OnlyPH.ll, f.Start.Val=f.Start.Val, VE.Start.Val=VE.Start.Val, delta.Start.Val=delta.Start.Val)))

      # Combine these results with those from previous iterations:
      VT.AnyAM.Result    = rbind(VT.AnyAM.Result,    VT.AnyAM.Result1)
      VT.OnlyPH.Result   = rbind(VT.OnlyPH.Result,   VT.OnlyPH.Result1)
      NoVT.AnyAM.Result  = rbind(NoVT.AnyAM.Result,  NoVT.AnyAM.Result1)
      NoVT.OnlyPH.Result = rbind(NoVT.OnlyPH.Result, NoVT.OnlyPH.Result1)
	  
	  i = i + 1 # Change the iteration number
    }
  }
}

###################################################################################################################################
# STEP 4: Choose the best iterations for each model, and compare them to determine which assumptions are most likely true/ untrue #
###################################################################################################################################

# Choose the iteration that best fit the data (i.e., with the lowest log likelihood):
    (Best.VT.AnyAM.Iteration    = VT.AnyAM.Result[which(VT.AnyAM.Result$ll       == min(VT.AnyAM.Result$ll,    na.rm=T)), 'Iteration'])
    (Best.VT.OnlyPH.Iteration   = VT.OnlyPH.Result[which(VT.OnlyPH.Result$ll     == min(VT.OnlyPH.Result$ll,   na.rm=T)), 'Iteration'])
    (Best.NoVT.AnyAM.Iteration  = NoVT.AnyAM.Result[which(NoVT.AnyAM.Result$ll   == min(NoVT.AnyAM.Result$ll,  na.rm=T)), 'Iteration'])
    (Best.NoVT.OnlyPH.Iteration = NoVT.OnlyPH.Result[which(NoVT.OnlyPH.Result$ll == min(NoVT.OnlyPH.Result$ll, na.rm=T)), 'Iteration'])

# Collect those iterations:
    Best.VT.AnyAM.Model    <- eval(parse(text=paste0("VT.AnyAM.It",    Best.VT.AnyAM.Iteration,    ".Fit")))
    Best.VT.OnlyPH.Model   <- eval(parse(text=paste0("VT.OnlyPH.It",   Best.VT.OnlyPH.Iteration,   ".Fit")))
    Best.NoVT.AnyAM.Model  <- eval(parse(text=paste0("NoVT.AnyAM.It",  Best.NoVT.AnyAM.Iteration,  ".Fit")))
    Best.NoVT.OnlyPH.Model <- eval(parse(text=paste0("NoVT.OnlyPH.It", Best.NoVT.OnlyPH.Iteration, ".Fit")))

# Compare the models to one another:
	mxCompare(Best.VT.AnyAM.Model, Best.NoVT.AnyAM.Model) # Tells us whether the trait is being influenced by parental effects
    mxCompare(Best.VT.AnyAM.Model, Best.VT.OnlyPH.Model)  # Tells us whether or not AM is being driven by phenotypic similarity

# Get the results from the best model:
	Best_Model <- Best.VT.AnyAM.Model # Which model performed best?

    Final_Result = data.frame(Model = Best_Model$name,
        VY = Best_Model$VY$values[1,1],
        VA = Best_Model$VY$values[1,1] - Best_Model$VF$values[1,1] - Best_Model$VE$values[1,1],
        VF = Best_Model$VF$values[1,1],
        VE = Best_Model$VE$values[1,1],
		delta = Best_Model$delta$values[1,1],
		Omega = Best_Model$Omega$values[1,1],
		k  = Best_Model$k$values[1,1],
        w  = Best_Model$w$result[1,1],
        f  = Best_Model$f$values[1,1],
		mu = Best_Model$mu$values[1,1],
        g  = Best_Model$g$values[1,1],
		Omega2mu = (Best_Model$Omega$values[1,1])^2 * Best_Model$mu$values[1,1],
        h2 = (Best_Model$VY$values[1,1] - Best_Model$VF$values[1,1] - Best_Model$VE$values[1,1])/Best_Model$VY$values[1,1],
        ll = Best_Model$output$Minus2LogLikelihood)

	# Final Result!
	Final_Result

	# If you prefer rounded numbers:
	(Final_Result_Rounded <- cbind(Final_Result$Model, round(Final_Result[,2:ncol(Final_Result)],4)))

	# If you prefer the variance components to be scaled such that they are out of 1:
	Final_Result[,c('VY','VA','VF','VE')] <- Final_Result[,c('VY','VA','VF','VE')] / Final_Result[,'VY']; Final_Result





