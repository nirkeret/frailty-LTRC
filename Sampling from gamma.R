
############################ Data Generation ####################################

n <- 5000 #desired sample size
nn <- 5000*5 #a bigger sample size out of which n observation will be randomly sampled. 

seed <- 1
set.seed(seed)


r <- 0.61 #type1 ceonsoring - observations with a greater observed time (event or censoring) will be censored at r.
g12 <- c(2,0.2,0.05)  #regression coefficients
g13 <- c(0.05,1)
g23 <- c(1,0.5)
theta <- 1  #frailty parameter
c12_1 <- 0.005 #constant hazard12 below recL
c12_2 <- 1 #constant hazard12 between recL and recU
c12_3 <- 1 #constant hazard12 above recU
c13_1 <- 0.5 #constant hazard13 below recL
c13_2 <- 1 #constant hazard13 between recL and recU
c13_3 <- 2 #constant hazard13 above recU
c23_1 <- 0 #constant hazard13 below lower23
c23_2 <- 1 #constant hazard13 between recL and recU
c23_3 <- 1 #constant hazard13 above recU

recL <- 0.05
lower23 <- 0.12
recU <- 0.15

#the "real" cummulative hazard functions
H012_real <- function(t) 
{
  ifelse(t < recL,c12_1*t,ifelse(t < recU,c12_1*recL + c12_2*(t - recL),c12_1*recL + c12_2*(recU - recL) + c12_3*(t - recU)))
}
H013_real <- function(t)
{
  ifelse(t < recL,c13_1*t,ifelse(t < recU,c13_1*recL + c13_2*(t - recL),c13_1*recL + c13_2*(recU - recL) + c13_3*(t - recU)))
}
H023_real <- function(t)
{
  ifelse(t < lower23,c23_1*t,ifelse(t < recU,c23_1*lower23 + c23_2*(t - lower23),c23_1*lower23 + c23_2*(recU - lower23) + c23_3*(t - recU)))
}  

######## creating reference data - mimicking the marginal death in the population below recL - 

nnn <- 100000
z1 <- runif(nnn) ; z2 <- runif(nnn) ; z3 <- runif(nnn)
z <- cbind(z1,z2,z3)

omega_ref <- rgamma(nnn,shape = 1/theta,scale = theta)
u13_ref <- runif(nnn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12_ref <- exp(Z12 %*% g12)
egz13_ref <- exp(Z13 %*% g13)

q1_ref <- egz12_ref*c12_1 + egz13_ref*c13_1
q2_ref <- egz12_ref*c12_2 + egz13_ref*c13_2
q3_ref <- egz12_ref*c12_3 + egz13_ref*c13_3

A13_recL_ref <- (c13_1*egz13_ref/(theta*q1_ref) * (exp(theta * q1_ref * recL) - 1))
S13_recL_ref <- exp(-omega_ref * A13_recL_ref)
A13_recU_ref <- A13_recL_ref + c13_2*egz13_ref/(theta*q2_ref)*exp(theta*q1_ref*recL)*(exp(theta*q2_ref*(recU-recL)) - 1)
S13_recU_ref <- exp(-omega_ref * A13_recU_ref)
s13_1_ref <- log(1 - (theta*q1_ref*log(u13_ref))/(omega_ref*c13_1*egz13_ref))/(theta*q1_ref)
s13_2_ref <- log(1 - (log(u13_ref)/omega_ref + A13_recL_ref)/(c13_2 * egz13_ref *exp(theta*q1_ref*recL)/(theta*q2_ref))) / (theta*q2_ref) + recL
s13_3_ref <- log(1 - (log(u13_ref)/omega_ref + A13_recU_ref)/(c13_3 * egz13_ref *exp(theta*(q1_ref*recL + q2_ref*(recU-recL)))/(theta*q3_ref))) / (theta*q3_ref) + recU

T13_ref <- ifelse(u13_ref > S13_recL_ref,s13_1_ref,ifelse(u13_ref > S13_recU_ref, s13_2_ref, s13_3_ref))
F13 <- ecdf(T13_ref)
H13 <- function(x) {-log(1 - F13(x))}
h13_approximate <- function(x,D) {(H13(x + D) - H13(x))/D}

n_grid <- 50
D <- 0.01
h13_grid_times <- seq(0,recL,length.out = n_grid)

h13_rates <- h13_approximate(h13_grid_times,D)


####main dataset sampling begins####
z1 <- runif(nn) ; z2 <- runif(nn) ; z3 <- runif(nn); z4 <- runif(nn)
z <- cbind(z1,z2,z3,z4)

omega <- rgamma(nn,shape = 1/theta,scale = theta)
u12 <- runif(nn)
u13 <- runif(nn)
u23 <- runif(nn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12 <- exp(Z12 %*% g12)
egz13 <- exp(Z13 %*% g13)

q1 <- egz12*c12_1 + egz13*c13_1
q2 <- egz12*c12_2 + egz13*c13_2
q3 <- egz12*c12_3 + egz13*c13_3


A12_recL <- (c12_1*egz12/(theta*q1) * (exp(theta * q1 * recL) - 1))
S12_recL <- exp(-omega * A12_recL)
A12_recU <- A12_recL + c12_2*egz12/(theta*q2)*exp(theta*q1*recL)*(exp(theta*q2*(recU-recL)) - 1)
S12_recU <- exp(-omega * A12_recU)
s12_1 <- log(1 - (theta*q1*log(u12))/(omega*c12_1*egz12))/(theta*q1)
s12_2 <- log(1 - (log(u12)/omega + A12_recL)/(c12_2 * egz12 *exp(theta*q1*recL)/(theta*q2))) / (theta*q2) + recL
s12_3 <- log(1 - (log(u12)/omega + A12_recU)/(c12_3 * egz12 *exp(theta*(q1*recL + q2*(recU-recL)))/(theta*q3))) / (theta*q3) + recU

T12 <- ifelse(u12 > S12_recL,s12_1,ifelse(u12 > S12_recU, s12_2, s12_3)) #healthy to diseased times

A13_recL <- (c13_1*egz13/(theta*q1) * (exp(theta * q1 * recL) - 1))
S13_recL <- exp(-omega * A13_recL)
A13_recU <- A13_recL + c13_2*egz13/(theta*q2)*exp(theta*q1*recL)*(exp(theta*q2*(recU-recL)) - 1)
S13_recU <- exp(-omega * A13_recU)
s13_1 <- log(1 - (theta*q1*log(u13))/(omega*c13_1*egz13))/(theta*q1)
s13_2 <- log(1 - (log(u13)/omega + A13_recL)/(c13_2 * egz13 *exp(theta*q1*recL)/(theta*q2))) / (theta*q2) + recL
s13_3 <- log(1 - (log(u13)/omega + A13_recU)/(c13_3 * egz13 *exp(theta*(q1*recL + q2*(recU-recL)))/(theta*q3))) / (theta*q3) + recU

T13 <- ifelse(u13 > S13_recL,s13_1,ifelse(u13 > S13_recU, s13_2, s13_3)) #death times - healthy observations

Z23 <- z[,c(1,4)]
egz23 <- exp(Z23 %*% g23)

A23_real <- function(t) {(exp(theta/(1+theta) * egz23 * H023_real(t)) - 1)/theta}
A23_lower <- A23_real(lower23)

S23_lower <- exp(-omega * A23_lower)
A23_recU <- A23_real(recU)
S23_recU <- exp(-omega * A23_recU)
S23_T12 <- exp(-omega * A23_real(T12))

Q <- (1+theta)/(theta*egz23) * log(theta*(A23_real(T12) - log(u23)/omega) + 1)
s23_1 <- Q/c23_1
s23_2 <- Q/c23_2 + (c23_2 - c23_1)*lower23/c23_2
s23_3 <- Q/c23_3 + (c23_2 - c23_1)*lower23/c23_3 + (c23_3-c23_2)*recU/c23_3

T23 <- ifelse(u23 > S23_lower,s23_1,ifelse(u23 > S23_recU, s23_2, s23_3))  #death times - diseased observations

R <- runif(nn,recL,recU) #Recruitment times


inout <- ifelse(T12 < T13, R < T23, R < T13) # observed in the sample
n.obs <- sum(inout)
Z12 <- Z12[inout,]; Z13 <- Z13[inout,]; Z23 <- Z23[inout,]   #retaining only the observed
T12 <- T12[inout]; T13 <- T13[inout]; T23 <- T23[inout]; R <- R[inout]

sample.indx <- sample(1:n.obs,n,replace=FALSE)  #sampling n observations out of the observed
Z12 <- Z12[sample.indx,]; Z13 <- Z13[sample.indx,]; Z23 <- Z23[sample.indx,]
T12 <- T12[sample.indx]; T13 <- T13[sample.indx]; T23 <- T23[sample.indx]; R <- R[sample.indx]


C <- rexp(n,2)  #independent censoring times

V <- pmin(T12,T13,R+C,r)  #first observed time
VpostR <- (V >= R)  #which observations are incident
delta1 <- (V == T12)  
delta1postR <- (V == T12) & VpostR
delta2 <- V == T13
W <- ifelse(!delta1,0,pmin(T23,R+C,r))  #second observed time 
delta3 <- as.vector((W == T23) & delta1)
