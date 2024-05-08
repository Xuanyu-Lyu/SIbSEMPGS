# The script is designed for using iterative math to check the algebraic expectations of the covariances

load("results100.rdata")

cor(results[[1]])
library(data.table)  

# Assuming results is your list of dataframes
cov_list <- lapply(results, cov, use="pairwise.complete.obs")

# Create an empty matrix to store the sums of the covariance matrices
sum_cov <- matrix(0, nrow = ncol(results[[1]]), ncol = ncol(results[[1]]))

# Sum up all the covariance matrices
for(cov_mat in cov_list) {
  sum_cov <- sum_cov + cov_mat
}

# Calculate the element-wise mean
mean_cov <- sum_cov / length(results)

#Example_Data  <- fread("Example_Data2.txt", header = T)

delta = .7
f = .2
MATE.COR = .3
k=.5
ObservedCovs = mean_cov

g = mean(c(ObservedCovs[5,c(6,7,8,11,12)],
           ObservedCovs[6,c(7,8,11,12)],
           ObservedCovs[7,c(8,9,10)],
           ObservedCovs[8,c(9,10)],
           ObservedCovs[9,c(10,11,12)],
           ObservedCovs[10,c(11,12)],
           ObservedCovs[11,c(12)])) 
omega = mean(c(ObservedCovs[3,c(5,6,9,10)],
               ObservedCovs[4,c(7,8,11,12)]))
w = (omega - delta*k - 2*delta*g)*2

covY1Nt2_expect = 2*delta*g+.5*delta*k+.5*w
covY1Nt2_obs = mean(c(ObservedCovs[1,9:12],
                      ObservedCovs[2,5:8]))

VY = mean(c(ObservedCovs[1,1],
            ObservedCovs[2,2],
            ObservedCovs[3,3],
            ObservedCovs[4,4]))
mu = ObservedCovs[3,4]/VY^2
VF = 2 * f^2 * VY * (1 + VY * mu)
covY1Y2_expect = 2*delta*w + 4*(delta^2)*g + delta^2*k + VF
covY1Y2_observed = ObservedCovs[1,2]


#save the observed covariance table into a csv file, round to 3 decimal places
write.csv(round(mean_cov,3), "ObservedCovs.csv", row.names = TRUE)

