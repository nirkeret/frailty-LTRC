library(statmod)

n <- 5000 #desired sample size
nn <- 5000*5 #a bigger sample size out of which n observation will be randomly sampled. 

r <- 0.61 #type1 ceonsoring - observations with a greater observed time (event or censoring) will be censored at r.
g12 <- c(2,0.2,0.05)  #regression coefficients
g13 <- c(0.05,1)
g23 <- c(1,0.5)
theta <- 1  #frailty parameter. for inverse Gaussian: 1/theta = variance. 
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


############sampling from the positive stable frailty ################
a.fun <- function(theta,alpha)
{
  num1<-sin( (1-alpha)*theta )
  num2<-sin( alpha*theta )^( alpha/(1-alpha) )
  den<-sin( theta)^( 1/(1-alpha) )
  out<-num1*num2/den
  out
}
ps.gen <- function(nobs,alpha)
{
  w<-rexp(nobs)
  theta<-runif(nobs)*pi
  out<-( a.fun(theta,alpha)/w )^( (1-alpha)/alpha )
  out
}
frailty.ps <- function(n, tau)
{
  alpha <- 1 - tau
  if(alpha > 0)
    omega <- ps.gen(n, alpha)
  if(alpha == 0)
    omega <- rep(1, n)
  return(omega)
}
## creating reference data
nnn <- 100000
z1 <- runif(nnn) ; z2 <- runif(nnn) ; z3 <- runif(nnn);
z <- cbind(z1,z2,z3)

omega_ref <- frailty.ps(nnn,tau)
u13_ref <- runif(nnn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12_ref <- exp(Z12 %*% g12)
egz13_ref <- exp(Z13 %*% g13)

q1_ref <- egz12_ref*c12_1 + egz13_ref*c13_1
q2_ref <- egz12_ref*c12_2 + egz13_ref*c13_2
q3_ref <- egz12_ref*c12_3 + egz13_ref*c13_3

A12_recL_ref <- c12_1*egz12_ref*q1_ref^((1-theta)/theta)*recL^(1/theta)
S12_recL_ref <- exp(-omega_ref * A12_recL_ref)
A12_recU_ref <- A12_recL_ref + c12_2*egz12_ref/q2_ref * ((q1_ref*recL + q2_ref*(recU-recL))^(1/theta) - (q1_ref*recL)^(1/theta))
S12_recU_ref <- exp(-omega_ref * A12_recU_ref)

tmp_ref <- q1_ref*recL + q2_ref*(recU-recL)
s12_1_ref <- (-log(u12_ref)/(c12_1*omega_ref))^theta * egz12_ref^(-theta) * q1_ref^(theta-1)
s12_2_ref <- ((q1_ref*recL)^(1/theta) - q2_ref/(c12_2*egz12_ref) * (log(u12_ref)/omega_ref + A12_recL_ref))^theta/q2_ref +recL*(1 - q1_ref/q2_ref)
s12_3_ref <- (tmp_ref^(1/theta) - q3_ref/(c12_3*egz12_ref) * (log(u12_ref)/omega_ref + A12_recU_ref))^theta/q3_ref - tmp_ref/q3_ref + recU

T12_ref <- ifelse(u12_ref > S12_recL_ref,s12_1_ref,ifelse(u12_ref > S12_recU_ref, s12_2_ref, s12_3_ref))

A13_recL_ref <- c13_1*egz13_ref*q1_ref^((1-theta)/theta)*recL^(1/theta)
S13_recL_ref <- exp(-omega_ref * A13_recL_ref)
A13_recU_ref <- A13_recL_ref + c13_2*egz13_ref/q2_ref * ((q1_ref*recL + q2_ref*(recU-recL))^(1/theta) - (q1_ref*recL)^(1/theta))
S13_recU_ref <- exp(-omega_ref * A13_recU_ref)

s13_1_ref <- (-log(u13_ref)/(c13_1*omega_ref))^theta * egz13_ref^(-theta) * q1_ref^(theta-1)
s13_2_ref <- ((q1_ref*recL)^(1/theta) - q2_ref/(c13_2*egz13_ref) * (log(u13_ref)/omega_ref + A13_recL_ref))^theta/q2_ref +recL*(1 - q1_ref/q2_ref)
s13_3_ref <- (tmp_ref^(1/theta) - q3_ref/(c13_3*egz13_ref) * (log(u13_ref)/omega_ref + A13_recU_ref))^theta/q3_ref - tmp_ref/q3_ref + recU

T13_ref <- ifelse(u13_ref > S13_recL_ref,s13_1_ref,ifelse(u13_ref > S13_recU_ref, s13_2_ref, s13_3_ref))

F13 <- ecdf(T13_ref)
H13 <- function(x) {-log(1 - F13(x))}
h13_apr <- function(x,D) {(H13(x + D) - H13(x))/D}

n_grid <- 50
D <- 0.01
h13_grid_times <- seq(0,recL,length.out = n_grid)

h13_grid <- h13_apr(h13_grid_times,D)

##creating the "main" sample:
z1 <- runif(nn) ; z2 <- runif(nn) ; z3 <- runif(nn); z4 <- runif(nn)
z <- cbind(z1,z2,z3,z4)

omega <- frailty.ps(nn,tau)
u12 <- runif(nn)
u13 <- runif(nn)
u23 <- runif(nn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12 <- exp(Z12 %*% g12)
egz13 <- exp(Z13 %*% g13)

q1 <- egz12*c12_1 + egz13*c13_1
q2 <- egz12*c12_2 + egz13*c13_2
q3 <- egz12*c12_3 + egz13*c13_3


A12_recL <- c12_1*egz12*q1^((1-theta)/theta)*recL^(1/theta)
S12_recL <- exp(-omega * A12_recL)
A12_recU <- A12_recL + c12_2*egz12/q2 * ((q1*recL + q2*(recU-recL))^(1/theta) - (q1*recL)^(1/theta))
S12_recU <- exp(-omega * A12_recU)

tmp <- q1*recL + q2*(recU-recL)
s12_1 <- (-log(u12)/(c12_1*omega))^theta * egz12^(-theta) * q1^(theta-1)
s12_2 <- ((q1*recL)^(1/theta) - q2/(c12_2*egz12) * (log(u12)/omega + A12_recL))^theta/q2 +recL*(1 - q1/q2)
s12_3 <- (tmp^(1/theta) - q3/(c12_3*egz12) * (log(u12)/omega + A12_recU))^theta/q3 - tmp/q3 + recU

T12 <- ifelse(u12 > S12_recL,s12_1,ifelse(u12 > S12_recU, s12_2, s12_3))

A13_recL <- c13_1*egz13*q1^((1-theta)/theta)*recL^(1/theta)
S13_recL <- exp(-omega * A13_recL)
A13_recU <- A13_recL + c13_2*egz13/q2 * ((q1*recL + q2*(recU-recL))^(1/theta) - (q1*recL)^(1/theta))
S13_recU <- exp(-omega * A13_recU)

s13_1 <- (-log(u13)/(c13_1*omega))^theta * egz13^(-theta) * q1^(theta-1)
s13_2 <- ((q1*recL)^(1/theta) - q2/(c13_2*egz13) * (log(u13)/omega + A13_recL))^theta/q2 +recL*(1 - q1/q2)
s13_3 <- (tmp^(1/theta) - q3/(c13_3*egz13) * (log(u13)/omega + A13_recU))^theta/q3 - tmp/q3 + recU

T13 <- ifelse(u13 > S13_recL,s13_1,ifelse(u13 > S13_recU, s13_2, s13_3))

Z23 <- z[,c(1,4)]
egz23 <- exp(Z23 %*% g23)

laplace_deriv <- function(x,theta) {-theta * exp(-x^theta) * x^(theta-1)}

A23_func_tosearch <- function(x,theta,egz,timepoint)
{
  tmp <- laplace_deriv(x,theta=theta) + exp(-H023_real(timepoint)*egz)
  tmp
}

A23_T12 <- rep(NA,nn)

for(i in 1:nn)
{
  A23_T12[i] <- uniroot(f=A23_func_tosearch,theta=theta,timepoint=T12[i],egz=egz23[i],interval = c(0,10000000),tol=10^(-20))$root
}

Q <- -log(-laplace_deriv(A23_T12 -log(u23)/omega,theta = theta)) / egz23 
H023_recL <- H023_real(recL)
H023_recU <- H023_real(recU)

s23_1 <- Q/c23_1
s23_2 <- Q/c23_2 + (c23_2 - c23_1)*lower23/c23_2
s23_3 <- Q/c23_3 + (c23_2 - c23_1)*lower23/c23_3 + (c23_3-c23_2)*recU/c23_3

T23 <- ifelse(Q < H023_recL, s23_1,ifelse(Q < H023_recU, s23_2, s23_3))

R <- runif(nn,recL,recU)

inout <- ifelse(T12 < T13, R < T23, R < T13) #observed in the sample
n.obs <- sum(inout)
Z12 <- Z12[inout,]; Z13 <- Z13[inout,]; Z23 <- Z23[inout,]
T12 <- T12[inout]; T13 <- T13[inout]; T23 <- T23[inout]; R <- R[inout]

sample.indx <- sample(1:n.obs,n,replace=FALSE)
Z12 <- Z12[sample.indx,]; Z13 <- Z13[sample.indx,]; Z23 <- Z23[sample.indx,]
T12 <- T12[sample.indx]; T13 <- T13[sample.indx]; T23 <- T23[sample.indx]; R <- R[sample.indx]


C <- rexp(n,2)

V <- pmin(T12,T13,R+C,r)
VpostR <- (V >= R)
delta1 <- (V == T12)  
delta1postR <- (V == T12) & VpostR
delta2 <- V == T13
W <- ifelse(!delta1,0,pmin(T23,R+C,r))
delta3 <- as.vector((W == T23) & delta1)


############sampling from the inverse Gaussian frailty ##############

quad_eq_solver <- function(a,b,c) {suppressWarnings((-b + sqrt(b^2 - 4*a*c))/(2*a))}


## creating a reference sample:
nnn <- 100000
z1 <- runif(nnn) ; z2 <- runif(nnn) ; z3 <- runif(nnn);
z <- cbind(z1,z2,z3)

omega_ref <- rinvgauss(nnn,mean = 1,shape = theta)
u12_ref <- runif(nnn) ; u13_ref <- runif(nnn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12_ref <- exp(Z12 %*% g12)
egz13_ref <- exp(Z13 %*% g13)

q1_ref <- egz12_ref*c12_1 + egz13_ref*c13_1
q2_ref <- egz12_ref*c12_2 + egz13_ref*c13_2
q3_ref <- egz12_ref*c12_3 + egz13_ref*c13_3

A13_recL_ref <- egz13_ref*c13_1*(q1_ref/(2*theta)*recL^2 + recL)
S13_recL_ref <- exp(-omega_ref * A13_recL_ref)
A13_recU_ref <- A13_recL_ref + egz13_ref*c13_2*((recL/theta * q1_ref + 1)*(recU-recL) + q2_ref/(2*theta)*(recU-recL)^2)
S13_recU_ref <- exp(-omega_ref * A13_recU_ref)

coef_13_1_a_ref <- egz13_ref*c13_1*q1_ref/(2*theta)
coef_13_1_b_ref <- egz13_ref*c13_1
coef_13_1_c_ref <- log(u13_ref)/omega_ref

coef_13_2_a_ref <- egz13_ref*c13_2*q2_ref/(2*theta)
coef_13_2_b_ref <- egz13_ref*c13_2*(recL/theta*q1_ref + 1)
coef_13_2_c_ref <- A13_recL_ref + log(u13_ref)/omega_ref

coef_13_3_a_ref <- egz13_ref*c13_3*q3_ref/(2*theta)
coef_13_3_b_ref <- egz13_ref*c13_3*(recL/theta*q1_ref + (recU-recL)/theta*q2_ref + 1)
coef_13_3_c_ref <- A13_recU_ref + log(u13_ref)/omega_ref

s13_1_ref <- quad_eq_solver(coef_13_1_a_ref,coef_13_1_b_ref,coef_13_1_c_ref)
s13_2_ref <- quad_eq_solver(coef_13_2_a_ref,coef_13_2_b_ref,coef_13_2_c_ref) + recL
s13_3_ref <- quad_eq_solver(coef_13_3_a_ref,coef_13_3_b_ref,coef_13_3_c_ref) + recU

T13_ref <- ifelse(u13_ref > S13_recL_ref,s13_1_ref,ifelse(u13_ref > S13_recU_ref, s13_2_ref, s13_3_ref))

F13 <- ecdf(T13_ref)
H13 <- function(x) {-log(1 - F13(x))}
h13_apr <- function(x,D) {(H13(x + D) - H13(x))/D}

n_grid <- 50
D <- 0.01
h13_grid_times <- seq(0,recL,length.out = n_grid)

h13_grid <- h13_apr(h13_grid_times,D)

## creating the "main" sample:

z1 <- runif(nn) ; z2 <- runif(nn) ; z3 <- runif(nn); z4 <- runif(nn)
z <- cbind(z1,z2,z3,z4)

omega <- rinvgauss(nn,mean = 1,shape = theta)
u12 <- runif(nn)
u13 <- runif(nn)
u23 <- runif(nn)

Z12 <- z[,1:3] ; Z13 <- z[,1:2]
egz12 <- exp(Z12 %*% g12)
egz13 <- exp(Z13 %*% g13)

q1 <- egz12*c12_1 + egz13*c13_1
q2 <- egz12*c12_2 + egz13*c13_2
q3 <- egz12*c12_3 + egz13*c13_3


A12_recL <- egz12*c12_1*(q1/(2*theta)*recL^2 + recL)
S12_recL <- exp(-omega * A12_recL)
A12_recU <- A12_recL + egz12*c12_2*((recL/theta * q1 + 1)*(recU-recL) + q2/(2*theta)*(recU-recL)^2)
S12_recU <- exp(-omega * A12_recU)

coef_12_1_a <- egz12*c12_1*q1/(2*theta)
coef_12_1_b <- egz12*c12_1
coef_12_1_c <- log(u12)/omega

coef_12_2_a <- egz12*c12_2*q2/(2*theta)
coef_12_2_b <- egz12*c12_2*(recL/theta*q1 + 1)
coef_12_2_c <- A12_recL + log(u12)/omega

coef_12_3_a <- egz12*c12_3*q3/(2*theta)
coef_12_3_b <- egz12*c12_3*(recL/theta*q1 + (recU-recL)/theta*q2 + 1)
coef_12_3_c <- A12_recU + log(u12)/omega

s12_1 <- quad_eq_solver(coef_12_1_a,coef_12_1_b,coef_12_1_c)
s12_2 <- quad_eq_solver(coef_12_2_a,coef_12_2_b,coef_12_2_c) + recL
s12_3 <- quad_eq_solver(coef_12_3_a,coef_12_3_b,coef_12_3_c) + recU


T12 <- ifelse(u12 > S12_recL,s12_1,ifelse(u12 > S12_recU, s12_2, s12_3))

A13_recL <- egz13*c13_1*(q1/(2*theta)*recL^2 + recL)
S13_recL <- exp(-omega * A13_recL)
A13_recU <- A13_recL + egz13*c13_2*((recL/theta * q1 + 1)*(recU-recL) + q2/(2*theta)*(recU-recL)^2)
S13_recU <- exp(-omega * A13_recU)

coef_13_1_a <- egz13*c13_1*q1/(2*theta)
coef_13_1_b <- egz13*c13_1
coef_13_1_c <- log(u13)/omega

coef_13_2_a <- egz13*c13_2*q2/(2*theta)
coef_13_2_b <- egz13*c13_2*(recL/theta*q1 + 1)
coef_13_2_c <- A13_recL + log(u13)/omega

coef_13_3_a <- egz13*c13_3*q3/(2*theta)
coef_13_3_b <- egz13*c13_3*(recL/theta*q1 + (recU-recL)/theta*q2 + 1)
coef_13_3_c <- A13_recU + log(u13)/omega

s13_1 <- quad_eq_solver(coef_13_1_a,coef_13_1_b,coef_13_1_c)
s13_2 <- quad_eq_solver(coef_13_2_a,coef_13_2_b,coef_13_2_c) + recL
s13_3 <- quad_eq_solver(coef_13_3_a,coef_13_3_b,coef_13_3_c) + recU

T13 <- ifelse(u13 > S13_recL,s13_1,ifelse(u13 > S13_recU, s13_2, s13_3))

Z23 <- z[,c(1,4)]
egz23 <- exp(Z23 %*% g23)

laplace_deriv <- function(x,theta)
{
  tmp <- sqrt(1 + 2/theta*x)
  -exp(theta*(1-tmp))/tmp
}
A23_func_tosearch <- function(x,theta,egz,timepoint)
{
  tmp <- laplace_deriv(x,theta = theta) + exp(-H023_real(timepoint)*egz)
  tmp
}

A23_T12 <- rep(NA,nn)

for(i in 1:nn)
{
  A23_T12[i] <- uniroot(f=A23_func_tosearch,theta = theta,timepoint=T12[i],egz=egz23[i],interval = c(0,10000000),tol=10^(-20))$root
}

Q <- -log(-laplace_deriv(A23_T12 -log(u23)/omega,theta = theta)) / egz23 
H023_recL <- H023_real(recL)
H023_recU <- H023_real(recU)

s23_1 <- Q/c23_1
s23_2 <- Q/c23_2 + (c23_2 - c23_1)*lower23/c23_2
s23_3 <- Q/c23_3 + (c23_2 - c23_1)*lower23/c23_3 + (c23_3-c23_2)*recU/c23_3

T23 <- ifelse(Q < H023_recL, s23_1,ifelse(Q < H023_recU, s23_2, s23_3))

R <- runif(nn,recL,recU)
 
inout <- ifelse(T12 < T13, R < T23, R < T13) #observed in the sample
n.obs <- sum(inout)
Z12 <- Z12[inout,]; Z13 <- Z13[inout,]; Z23 <- Z23[inout,]
T12 <- T12[inout]; T13 <- T13[inout]; T23 <- T23[inout]; R <- R[inout]

sample.indx <- sample(1:n.obs,n,replace=FALSE)  #sampling n observations out of nn
Z12 <- Z12[sample.indx,]; Z13 <- Z13[sample.indx,]; Z23 <- Z23[sample.indx,]
T12 <- T12[sample.indx]; T13 <- T13[sample.indx]; T23 <- T23[sample.indx]; R <- R[sample.indx]


C <- rexp(n,2)

V <- pmin(T12,T13,R+C,r)
VpostR <- (V >= R)
delta1 <- (V == T12)  
delta1postR <- (V == T12) & VpostR
delta2 <- V == T13
W <- ifelse(!delta1,0,pmin(T23,R+C,r))
delta3 <- as.vector((W == T23) & delta1)
