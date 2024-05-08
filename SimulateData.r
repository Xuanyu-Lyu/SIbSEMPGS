source("SIMULATE_DAT_GEN.R")

num.cvs = 500

RUN.MARKERS <- FALSE #whether to only consider GRMs built from CVs (FALSE) or both CVs and SNPs (TRUE)

MIN.MAF <- .1; MAX.MAF <- .50#AM Simulation wildcards

MAF.VECTOR <- runif(num.cvs,MIN.MAF,MAX.MAF) #Can change the distribution of MAFs here

GENTP.VAR <- MAF.VECTOR*(1-MAF.VECTOR)*2

ALPHA.VECTOR <- sample(c(-1,1),num.cvs,replace=TRUE)*sqrt(1/(num.cvs*GENTP.VAR)) #Can change the distribution of effect sizes here - fixed f'n of MAF

CV.INFO <- data.frame(MAF=MAF.VECTOR,alpha=ALPHA.VECTOR)

data_list <- AM.SIMULATE(
    CV.INFO = CV.INFO, 
    H2.T0 = .49, 
    NUM.GENERATIONS = 20, 
    POP.SIZE = 20000, 
    MATE.COR = .3, 
    AVOID.INB = TRUE, 
    SAVE.EACH.GEN = TRUE, 
    SAVE.COVS = TRUE, 
    SEED = 62, 
    VF.T0 = .08, 
    PROP.H2.LATENT = 0, 
    Unequal_AM = TRUE)

data_df <- data_list$HISTORY$PHEN[[11]] |> as.data.frame()


all_dadID <- unique(data_df[,2])
all_trios <- data.frame(matrix(ncol = ncol(data_df), nrow = 0))
colnames(all_trios) <- colnames(data_df)
for (i in 1:length(all_dadID)){
    dadID <- all_dadID[i]
    trio_df <- data_df[data_df[,2] == dadID,] 
    if (nrow(trio_df) >= 2){
        all_trios <- rbind(all_trios,trio_df[sample.int(nrow(trio_df),2),])
    } else {
       next
    }
    
}


num_rows <- nrow(all_trios)/2

df_final <- data.frame(Yo1 = rep(NA, num_rows),
                       Yo2 = rep(NA, num_rows),
                       Yp = rep(NA, num_rows),
                       Ym = rep(NA, num_rows),
                       Tp1 = rep(NA, num_rows),
                       NTp1 = rep(NA, num_rows),
                       Tm1 = rep(NA, num_rows),
                       NTm1 = rep(NA, num_rows),
                       Tp2 = rep(NA, num_rows),
                       NTp2 = rep(NA, num_rows),
                       Tm2 = rep(NA, num_rows),
                       NTm2 = rep(NA, num_rows))
df_final$Yo1 <- all_trios$Y[seq(1, nrow(all_trios), 2)]
df_final$Yo2 <- all_trios$Y[seq(2, nrow(all_trios), 2)]
df_final$Yp <- all_trios$YP[seq(1, nrow(all_trios), 2)]
df_final$Ym <- all_trios$YM[seq(1, nrow(all_trios), 2)]
df_final$Tp1 <- all_trios$TPO[seq(1, nrow(all_trios), 2)]
df_final$NTp1 <- all_trios$NTPO[seq(1, nrow(all_trios), 2)]
df_final$Tm1 <- all_trios$TMO[seq(1, nrow(all_trios), 2)]
df_final$NTm1 <- all_trios$NTMO[seq(1, nrow(all_trios), 2)]
df_final$Tp2 <- all_trios$TPO[seq(2, nrow(all_trios), 2)]
df_final$NTp2 <- all_trios$NTPO[seq(2, nrow(all_trios), 2)]
df_final$Tm2 <- all_trios$TMO[seq(2, nrow(all_trios), 2)]
df_final$NTm2 <- all_trios$NTMO[seq(2, nrow(all_trios), 2)]

write.table( df_final,"Example_Data2.txt", row.names = FALSE)

# a loop to run the simulation 100 times and save the results in a list and write a rdata rile

results <- list()
for(k in 1:100){
    data_list <- AM.SIMULATE(
        CV.INFO = CV.INFO, 
        H2.T0 = .49, 
        NUM.GENERATIONS = 20, 
        POP.SIZE = 20000, 
        MATE.COR = .3, 
        AVOID.INB = TRUE, 
        SAVE.EACH.GEN = TRUE, 
        SAVE.COVS = TRUE, 
        SEED = k, 
        VF.T0 = .08, 
        PROP.H2.LATENT = 0, 
        Unequal_AM = TRUE)

    data_df <- data_list$HISTORY$PHEN[[11]] |> as.data.frame()
    
    all_dadID <- unique(data_df[,2])
    all_trios <- data.frame(matrix(ncol = ncol(data_df), nrow = 0))
    colnames(all_trios) <- colnames(data_df)
    for (i in 1:length(all_dadID)){
        dadID <- all_dadID[i]
        trio_df <- data_df[data_df[,2] == dadID,] 
        if (nrow(trio_df) >= 2){
            all_trios <- rbind(all_trios,trio_df[sample.int(nrow(trio_df),2),])
        } else {
        next
        }
        
    }
    num_rows <- nrow(all_trios)/2
        df_final <- data.frame(Yo1 = rep(NA, num_rows),
                            Yo2 = rep(NA, num_rows),
                            Yp = rep(NA, num_rows),
                            Ym = rep(NA, num_rows),
                            Tp1 = rep(NA, num_rows),
                            NTp1 = rep(NA, num_rows),
                            Tm1 = rep(NA, num_rows),
                            NTm1 = rep(NA, num_rows),
                            Tp2 = rep(NA, num_rows),
                            NTp2 = rep(NA, num_rows),
                            Tm2 = rep(NA, num_rows),
                            NTm2 = rep(NA, num_rows))
        df_final$Yo1 <- all_trios$Y[seq(1, nrow(all_trios), 2)]
        df_final$Yo2 <- all_trios$Y[seq(2, nrow(all_trios), 2)]
        df_final$Yp <- all_trios$YP[seq(1, nrow(all_trios), 2)]
        df_final$Ym <- all_trios$YM[seq(1, nrow(all_trios), 2)]
        df_final$Tp1 <- all_trios$TPO[seq(1, nrow(all_trios), 2)]
        df_final$NTp1 <- all_trios$NTPO[seq(1, nrow(all_trios), 2)]
        df_final$Tm1 <- all_trios$TMO[seq(1, nrow(all_trios), 2)]
        df_final$NTm1 <- all_trios$NTMO[seq(1, nrow(all_trios), 2)]
        df_final$Tp2 <- all_trios$TPO[seq(2, nrow(all_trios), 2)]
        df_final$NTp2 <- all_trios$NTPO[seq(2, nrow(all_trios), 2)]
        df_final$Tm2 <- all_trios$TMO[seq(2, nrow(all_trios), 2)]
        df_final$NTm2 <- all_trios$NTMO[seq(2, nrow(all_trios), 2)]
    results[[k]] <- df_final
    cat("Simulation ", k, " completed\n")

}
    # write the list into a rdata file
    save(results, file = "results.rdata")
