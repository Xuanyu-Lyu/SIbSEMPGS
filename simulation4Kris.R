# Step 1: simulate a group of data
source("GeneEvolve/SIMULATE_DAT_GEN.R")

# Initialize the genetic information
num.cvs = 100

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
AM <- c(0,.5)


for (i in 1:length(AM)){
    # a_true <- getAfromR2pgs(v_r2pgs[i], h2)
    # cat(a_true,"\n")
    # PLat <- getPLatfromR2pgs(v_r2pgs[i], h2)
    # run the simulation to get some datasets

    results <- list()
    for(k in 1:20){
        data_list <- AM.SIMULATE(
            CV.INFO = CV.INFO, 
            H2.T0 = .5, 
            NUM.GENERATIONS = 20, 
            POP.SIZE = 50000, 
            MATE.COR = AM[i], 
            AVOID.INB = TRUE, 
            SAVE.EACH.GEN = TRUE, 
            SAVE.COVS = TRUE, 
            SEED = k, 
            VF.T0 = 0, 
            PROP.H2.LATENT = 0, 
            Unequal_AM = FALSE)

        data_df <- data_list$HISTORY 
        results[[k]] <- data_df
        cat("AM =", AM[i], "Simulation ", k, " completed\n")
        
    }
    save(results, file = paste0("AM", AM[i],"results.rdata"))
}

# check the data
load("AM0results.rdata")
results[[1]] |> str()
offspring_example <-  results[[1]]["XO"]
offspring_example$XO[20] |> as.data.frame() |> head()
results[[1]]$PHEN[[20]] |> as.data.frame() |> psych::describe() |> print()
data_AM0 <- list()
for (i in 1:20){
    datai <- cbind(as.data.frame(results[[i]]$PHEN[[20]])$Y,as.data.frame(results[[i]]["XO"]$XO[20]))
    colnames(datai)[1] <- c("Y")
    data_AM0[[i]] <- datai
}
str(data_AM0)
data_AM0[[1]] |> head()

saveRDS(data_AM0, "data_AM0.rds")

load("AM0.5results.rdata")
data_AM0.5 <- list()
for (i in 1:20){
    datai <- cbind(as.data.frame(results[[i]]$PHEN[[20]])$Y,as.data.frame(results[[i]]["XO"]$XO[20]))
    colnames(datai)[1] <- c("Y")
    data_AM0.5[[i]] <- datai
}
str(data_AM0.5)
data_AM0.5[[1]] |> head()
saveRDS(data_AM0.5, "data_AM0.5.rds")


data_AM0.5 <- list()
for (i in 1:20){
    datai <- cbind(as.data.frame(results[[i]]$PHEN[[1]])$Y,as.data.frame(results[[i]]["XO"]$XO[1]))
    colnames(datai)[1] <- c("Y")
    data_AM0.5[[i]] <- datai
}
str(data_AM0.5)
data_AM0.5[[1]] |> var()
saveRDS(data_AM0.5, "data_AM0.5_firstG.rds")

