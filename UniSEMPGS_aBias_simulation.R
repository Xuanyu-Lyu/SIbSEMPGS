# This is a script to analyse how misspecified a can lead to bias in other parameters. 
# Author: Xuanyu Lyu
# Date: 08/24/2024

# meet with Dinka 08/27 :
# 1. level of true r2pgs: .025/.05/.1
# 2.1) level of bias r2pgs: .025- .01/.04/.2
#   2) level of bias r2pgs: .05- .025/.05/.075/.2
#   3) level of bias r2pgs: .1- .05/.15/.2
# 3. make plots: x-axis: level of bias r2pgs, y-axis: estimates of the targeted parameters
# 4. make plots: vary a from .5 to .8 in .1 increments and draft the plots
# Dinka's bd: 09/10



# Step 1: simulate a group of data
source("SIMULATE_DAT_GEN.R")

# Initialize the genetic information
num.cvs = 200

RUN.MARKERS <- FALSE #whether to only consider GRMs built from CVs (FALSE) or both CVs and SNPs (TRUE)

MIN.MAF <- .1; MAX.MAF <- .50#AM Simulation wildcards

MAF.VECTOR <- runif(num.cvs,MIN.MAF,MAX.MAF) #Can change the distribution of MAFs here

GENTP.VAR <- MAF.VECTOR*(1-MAF.VECTOR)*2

ALPHA.VECTOR <- sample(c(-1,1),num.cvs,replace=TRUE)*sqrt(1/(num.cvs*GENTP.VAR)) #Can change the distribution of effect sizes here - fixed f'n of MAF

CV.INFO <- data.frame(MAF=MAF.VECTOR,alpha=ALPHA.VECTOR)

getPLatfromR2pgs <- function(r2pgs,h2){
    return((h2-r2pgs)/h2)
}

getAfromR2pgs <- function(r2pgs, h2){
    return(sqrt(h2-r2pgs))
}
v_r2pgs <- c(.025,.05,.1)
h2 <- .49


# get all the simulated data
for (i in 1:length(v_r2pgs)){
    a_true <- getAfromR2pgs(v_r2pgs[i], h2)
    cat(a_true,"\n")
    PLat <- getPLatfromR2pgs(v_r2pgs[i], h2)
    # run the simulation to get some datasets

    results <- list()
    for(k in 1:100){
        data_list <- AM.SIMULATE(
            CV.INFO = CV.INFO, 
            H2.T0 = h2, 
            NUM.GENERATIONS = 20, 
            POP.SIZE = 20000, 
            MATE.COR = .3, 
            AVOID.INB = TRUE, 
            SAVE.EACH.GEN = TRUE, 
            SAVE.COVS = TRUE, 
            SEED = k, 
            VF.T0 = .08, 
            PROP.H2.LATENT = PLat, 
            Unequal_AM = FALSE)

        data_df <- data_list$HISTORY$PHEN[[20]] |> as.data.frame()
        results[[k]] <- data_df
        cat("r2pgs =", v_r2pgs[i], "Simulation ", k, " completed\n")
        
    }
    save(results, file = paste0("r2pgs", v_r2pgs[i],"results_biasA100.rdata"))
}




    a_true <- sqrt(.49*.7)
    # write the list into a rdata file
    







