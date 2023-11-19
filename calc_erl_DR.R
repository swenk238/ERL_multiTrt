#1. Use original Ts and X's to compute the Z_Ts and Z_Xs through ERL
#2. Will get latent level with MVN mean 0
#3. Get Z_T | Z_X -- propensity score Z_R1(Z_T, Z_X),on the latent level
#4. Fit H ~ RAW SCALED T and Z_R -- get beta_hat's -- $\mu = \beta_0 + \beta_1T_1  + \beta_3T_2 + \beta_5Z_{R^2} + \beta_6 T_1T_2Z_R$
#5. Compute Z_R2(Z_T^*) -- counterfactual Z_R1
#6. Then plug them all back into (4) to get the dose response using beta's estimated from (4) to get expected H's
#$\mu^* = \beta_0 + \beta_1T^*_1  + \beta_3T^*_2 + \beta_5Z_{Z_{R2}} + \beta_6 T^*_1T^*_2Z_{R2}$
#7. Then plot 3D plot of E(H | Z_R2, T1, T2) \mu^* against T^*_1 and T^*_2

library(dplyr)
Rcpp::sourceCpp('cppFunc.cpp')

wd <- "~/GitHub/ERL_multiTrt"
setwd(wd)

# Get LGPS, Z_R for all the countries
load("lacsh_scaledZ.RData")

## FROM calc_lacsh_scaledZ.R
## using scaled Z
mu_hat_mat2 = list()
Sig_hat_mat2 = list()
lgps_list2 = lapply(1:NSCAN, function(x)matrix(NA, nrow=n))

for(nscan in 1:NSCAN){
  mu_hat_mat2[[nscan]] = MatMult3(as.matrix(C_list[[nscan]][c(1,2),-c(1,2)]), solve(C_list[[nscan]][-c(1,2),-c(1,2)]), as.matrix(scaled_Z_list[[nscan]][-c(1,2),]))
  Sig_hat_mat2[[nscan]] = as.matrix(C_list[[nscan]][c(1,2),c(1,2)]) - MatMult3(as.matrix(C_list[[nscan]][c(1,2),-c(1,2)]), solve(C_list[[nscan]][-c(1,2),-c(1,2)]), (C_list[[nscan]][-c(1,2),c(1,2)]))
  
  # LGPS for two treatments Z_T1,T2 | Z_X
  for(i in 1:N){
    lgps_list2[[nscan]][i,] = dmvnorm(scaled_Z_list[[nscan]][c(1,2),i], mean = mu_hat_mat[[nscan]][,i], sigma= Sig_hat_mat[[nscan]])
  }
}

# extract median of H's, and raw scaled_T and do regression 
H_med <- readRDS("H_med_model_B1.RDS")

DR_dat <- cbind(H_med, TRT1, TRT2, lgps_list2[[1]]) %>% data.frame()
colnames(DR_dat) <- c('H_med', 'T1', 'T2', 'Z_R')

# getting estimates
fit.lm <- lm(H_med ~ T1 + T2 + Z_R + T1:T2:Z_R, data = DR_dat) 
fit.lm %>% summary()

fit.coef <- fit.lm$coefficients

## specified T's 
seqT1 = seq(range(TRT1)[1],range(TRT1)[2], length.out=30)
seqT2 = seq(range(TRT2)[1],range(TRT2)[2], length.out=30)

## Z_T1T2 -- extracting Z_t1 and Z_t2 from scaled_Z_list 
scaled_Z_list[[2]][1:2,]
Z_T1 = lapply(scaled_Z_list[1:10],'[',1,) %>% do.call(rbind,.) %>% apply(.,2,median)
Z_T2 = lapply(scaled_Z_list[1:10],'[',2,) %>% do.call(rbind,.) %>% apply(.,2,median)

rank(Z_T1)

Z_T1T2 = cbind(Z_T1, Z_T2)

## mapping specified Ts to Z_T's to get Z_specified T
Z_T1star = map_rank(Z_T1, seqT1) %>% as.matrix()
Z_T2star = map_rank(Z_T2, seqT2) %>% as.matrix()

# Propensity score for two treatments 
Z_r_dualtrt = matrix(NA,nrow=N, ncol=length(seqT1))

for(i in 1:N){
  Z_r_dualtrt[i,] = dmvnorm(cbind(Z_T1star, Z_T2star), mean = mu_hat_mat[[nscan]][,i], sigma= Sig_hat_mat[[nscan]])
}

len = length(Z_T1star)
Z_r_dualtrt = lapply(1:N, function(x)matrix(NA,ncol=len, nrow=len))
for(i in 1:2){ # should be 1:N
  for (k in 1:len){
    for (j in 1:len){
      Z_r_dualtrt[[i]][j,k] = dmvnorm(c(Z_T1star[k,1], Z_T2star[j,1]), mean = mu_hat_mat[[nscan]][,i], sigma= Sig_hat_mat[[nscan]])
    }
  }
}
# cbind(sort(Z_T1), sort(seqT1), map_rank(Z_T1, seqT1)) %>% View() # looks right

# test code to map the ranks ----------------------------------------------
x = c(-3, 0, 3) %>% sort() # latent T
y = c(-2,-1,1,4,6,2,7) %>% sort() # specified T

map_rank <- function(x,y){ # it'll correspond to the new values for SORTED original Y
  x <- sort(x)
  y <- sort(y)
  p = length(x)-1
  rank_y = matrix(NA, ncol=length(y), nrow=p)
  for(j in 1:p){
    rank_y[j,] = ifelse(y >= x[j] & y < x[j+1], rank(x)[j], # if y value is in between x1 and x2, take x1's rank
                        ifelse(y > x[j+1], rank(x)[j+1], # if y value is bigger than largest x then take x
                               ifelse(y < x[1], rank(x)[1], 0)))# if y is smaller than smallest x then take x
  }
  
  rank_y_vec <- apply(rank_y, 2,max)
  print(x[rank_y_vec])
}
map_rank(x,y) # answer should be c(-3,-3,0,0,3,3,3)




# from sbgcop.mcmc()
Y <- testZ_T
Y <- as.matrix(Y)
n <- dim(Y)[1]
p <- dim(Y)[2]
R <- NULL
for (j in 1:p) {
  R <- cbind(R, match(Y[, j], sort(unique(Y[, j]))))
}
Rlevels <- apply(R, 2, max, na.rm = TRUE)
Ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")
N <- apply(!is.na(Ranks), 2, sum)
U <- t(t(Ranks)/(N + 1))
Z <- qnorm(U)

n = nrow(Z_T1T2)
p = ncol(Z_T1T2)

########## constraints
## ranking the values of Y for each column
R<-NULL 
for(j in 1:p) {
  R<-cbind(R, match(Z_T1T2[,j],sort(unique(Z_T1T2[,j])))) 
}


for(r in sort(unique(R[,j]))){
  ir<- (1:n)[R[,j]==r & !is.na(R[,j])]
  lb<-suppressWarnings( max(Z[ R[,j]<r,j],na.rm=T) ) #extracting the Z rows where the values are less than the specified rank in R  
  ub<-suppressWarnings( min(Z[ R[,j]>r,j],na.rm=T) )
  Z[ir,j]<-qnorm(runif(length(ir), pnorm(lb,muj[ir],sdj),pnorm(ub,muj[ir],sdj)),muj[ir],sdj)
}
ir<-(1:n)[is.na(R[,j])]






# calculate dose response  ------------------------------------------------
## dose response at each specified trt level 
DR_list = lapply(1:N, function(x)matrix(NA,ncol=len, nrow=len)) 
for (i in 1:2){ #should be 1:N
  for (k in 1:len){
    for (j in 1:len){
    DR_list[[i]][j,k] = fit.coef%*%c(1,seqT1[j], seqT2[j], Z_r_dualtrt[[i]][j,k], seqT1[j]*seqT2[j]*Z_r_dualtrt[[i]][j,k])
  }
  }
}

DR_mean = matrix(NA, ncol=len, nrow=len)
for (j in 1:len){
  DR_mean[,j] = apply(sapply(DR_list[1:2],'[',,j),1,median)
}

save(DR_mean, seqT1, seqT2, file="ERL_dualtrt_DR.RData")
