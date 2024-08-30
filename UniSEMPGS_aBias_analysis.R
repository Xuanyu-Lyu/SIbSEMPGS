# Step 2: fit the model with different levels of a
source("UniSEMPGS_OpenMx_aBias.R")
# load the data
#load("results_biasA100.rdata")

# run a test trial to debug
# test_df = results[[1]][,c("NTMO","TMO","NTPO","TPO","YM","YP","Y")]
# colnames(test_df) <- c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")

# test_result <- fitDifferentA(test_df, a_true, max.cores = 2)
# test_result


# test done, now fit the model to all simulated datasets



# a function to average across all the result dfs
averageResult <- function(results){
    result_df <- results[[1]]
    for(k in 2:length(results)){
        result_df <- result_df + results[[k]]
    }
    result_df <- result_df/length(results)
    return(result_df)
}

getAfromR2pgs <- function(r2pgs, h2){
    return(sqrt(h2-r2pgs))
}

v_r2pgs <- c(.025,.05,.1)
v_a <- sapply(v_r2pgs, getAfromR2pgs, h2 = .49) 
# r2pgs = 0.025
# r2pgs_biased = .01/.04/.2
v_r2pgs_biased <- matrix(c(c(.025,.01,.04,.2),
                               c(.05,.025,.075,.2),
                               c(.1,.05,.15,.2)), ncol = 4, byrow = TRUE)
v_a_biased <- apply(v_r2pgs_biased, c(1,2), getAfromR2pgs, h2 = .49)
l_r2pgs_fitsum <- list()

v_a

for (i in 1:length(v_r2pgs)){
    l_r2pgs_fitsum[[as.character(v_r2pgs[i])]] <- list()
    load(paste0("r2pgs", v_r2pgs[i],"results_biasA100.rdata"))
    for(j in 1:length(v_r2pgs_biased[i,])){
        l_r2pgs_fitsum[[as.character(v_r2pgs[i])]][[as.character(v_r2pgs_biased[i,j])]] <- list()
        for(k in 1:100){
        data_df <- results[[k]]
        data_df <- data_df[,c("NTMO","TMO","NTPO","TPO","YM","YP","Y")]
        colnames(data_df) <- c("NTm","Tm","NTp","Tp","Ym","Yp","Yo")
        l_r2pgs_fitsum[[as.character(v_r2pgs[i])]][[as.character(v_r2pgs_biased[i,j])]][[k]] <- fitDifferentA(data_df, a = v_a_biased[i,j] , max.cores = 2)[,-1]
        cat("r2pgs =", v_r2pgs[i],"biased r2pgs =",v_a_biased[i,j], "data", k, " completed\n")
        }
    }
   
}
save(l_r2pgs_fitsum, file = "l_r2pgs_fitsum.rdata")

# Analyze the results
# load the results
load("l_r2pgs_fitsum.rdata")
# extract the desired parameters from the 100 fits and put them in a vector
extractParam <- function(List, param){
    v_param <- c()
    for(i in 1:length(List)){
        v_param[i] <- List[[i]][1,param]
    }
    return(v_param)
}


# true values
#Omega 0.3210373
#Gamma 0.4903925
#VY 1.723142
#g 0.01794366
#h 0.04186854
#w 0.1669394
#v 0.2550041
#VF 0.1792068


# make figure for the VF parameter
# Create column names
col_names <- c(paste0("r2pgs_bias", round(v_r2pgs_biased[1,2],2)),
               paste0("r2pgs_bias", round(v_r2pgs_biased[1,3],2)),
               paste0("r2pgs_bias", round(v_r2pgs_biased[1,4],2)))

# Create the data frame
df_VF <- data.frame(
    a_bias1 = extractParam(l_r2pgs_fitsum[["0.025"]][["0.01"]], "VF"),
    a_bias2 = extractParam(l_r2pgs_fitsum[["0.025"]][["0.04"]], "VF"),
    a_bias3 = extractParam(l_r2pgs_fitsum[["0.025"]][["0.2"]], "VF")
)
colnames(df_VF) <- col_names
head(df_VF)
# true_values for biased analysis (replace with actual values)
true_values <- 0.1792068

# Convert df_VF to long format
df_long_biased <- tidyr::pivot_longer(df_VF, cols = everything(), names_to = "Variable", values_to = "Value")
df_long_biased$Index <- 1:nrow(df_long_biased)

# Calculate means for biased variables
means_biased <- aggregate(Value ~ Variable, data = df_long_biased, FUN = mean)
library(ggplot2)
# Create the plot
# Create the plot
g_biased <- ggplot(df_long_biased, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values), color = "red") +
  geom_hline(data = means_biased, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  #theme_minimal() +
  annotate("text", x = Inf, y = true_values, label = round(true_values, 3), 
           color = "red", hjust = 1.1, vjust = -0.5, size = 3) +
  # Add annotations for mean values
  geom_text(data = means_biased, aes(x = Inf, y = Value, label = round(Value, 3)), 
            color = "blue", hjust = 1.1, vjust = 1.5, size = 3)

# Display the plot
g_biased
ggsave("VF.png", g_biased, width = 8, height = 8, type = "cairo-png", dpi = 400)


# make figure for the VF parameter
# Create column names
col_names <- c(paste0("a_bias", round(v_a_biased[2,2],2)),
               paste0("a_bias", round(v_a_biased[2,3],2)),
               paste0("a_bias", round(v_a_biased[2,4],2)))

# Create the data frame
df_VF <- data.frame(
    a_bias1 = extractParam(l_r2pgs_fitsum[["0.05"]][["0.025"]], "VF"),
    a_bias2 = extractParam(l_r2pgs_fitsum[["0.05"]][["0.075"]], "VF"),
    a_bias3 = extractParam(l_r2pgs_fitsum[["0.05"]][["0.2"]], "VF")
)
colnames(df_VF) <- col_names
head(df_VF)
# true_values for biased analysis (replace with actual values)
true_values <- 0.1792068

# Convert df_VF to long format
df_long_biased <- tidyr::pivot_longer(df_VF, cols = everything(), names_to = "Variable", values_to = "Value")
df_long_biased$Index <- 1:nrow(df_long_biased)

# Calculate means for biased variables
means_biased <- aggregate(Value ~ Variable, data = df_long_biased, FUN = mean)
library(ggplot2)
# Create the plot
# Create the plot
g_biased <- ggplot(df_long_biased, aes(x = Index, y = Value)) +
  geom_point() +
  geom_hline(aes(yintercept = true_values), color = "red") +
  geom_hline(data = means_biased, aes(yintercept = Value), color = "blue", linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free") +
  #theme_minimal() +
  annotate("text", x = Inf, y = true_values, label = round(true_values, 3), 
           color = "red", hjust = 1.1, vjust = -0.5, size = 3) +
  # Add annotations for mean values
  geom_text(data = means_biased, aes(x = Inf, y = Value, label = round(Value, 3)), 
            color = "blue", hjust = 1.1, vjust = 1.5, size = 3)

# Display the plot
g_biased
ggsave("VF.png", g_biased, width = 8, height = 8, type = "cairo-png", dpi = 400)
