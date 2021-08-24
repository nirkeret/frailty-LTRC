# frailty-LTRC
The "Estimation Procedure" file contains the code for running the algorithm. It can be given real data or simulated data, such as by the means of the codes provided for sampling from the gamma or IG/PS frailties.
Before running the estimation procedure, the cpp file in this directory should be placed in the same working directory (or elsewhere, but the cpp file should be read before executing the estimation procedure).
Detailed descriptions for how to run the estimation code are provided as comments on the code itself, but for convenience it is repeated here:

## Estimation Function Arguments 

z12 - covariate matrix (n_subjects X n_covariates) for the healthy-disease process.

z13 - covariate matrix (n_subjects X n_covariates) for the healthy-dead process.

z23 - covariate matrix (n_subjects X n_covariates) for the disease-dead process.

obs_times1 - minimum of disease time, death time(from healthy state) or censoring time.

obs_times2 - minimum of death time (from diseased state) and censoring. If not applicable (if the person died or was censored directly from healthy state) should get 0.

R - recruitment(delated entry) age. Defaults to 0 (no delayed entry).

delta1 - TRUE if obs_times1 = disease time, FALSE otherwise.

delta2 - TRUE if obs_times1 = death time, FALSE otherwise.

delta3 - TRUE if obs_times2 = death time, FALSE otherwise.

RecL - minimal recruitment age. default value is 0.

marginal_13_haz - the marginal hazard rates, across a grid until recL. defaults to NULL. has to be provided if recL is not 0.

marginal_13_grid - the marginal hazard grid times corresponding to the marginal hazard rates. defaults to NULL. has to be provided if recL is not 0.

n_Bootstrap - number of required weighted bootstrap samples. Defaults to none. Weights are randomly sampled from the standard exponential distribution.

Tol - tolerance for the iterative algorithm. If |new_loglikelikehood - last_loglikelihood|/|last_loglikelihood| < tol, the algorithm stops and declares convergence. Defaults to 0.001.

max_Iter - maximal number of iterations, after which the algorithm stops and declares non-convergence. Defaults to 100.

initial_pars - a vector with starting points for theta+regression parameters. Defaults to naive Cox regression with LT. Theta value defaults to 0.1.

initial_H - a list with 3 starting functions for H012,H013,H023. Defaults to naive Breslow estimator with LT. 

seed - controlling the seed for the weights generation. Defaults to 1.

save_address - a local address to save interim results. After each bootstrap round the results will be saved to the address provided. The address should direct to a folder, where 5 files will be saved: H012,H013,H023,parameters,convergence. **The first row in each file corresponds to the original data point estimates. The following rows correspond to the point estimates based on bootstraped data.**

check_times - a vector of times to be inserted into the H0 functions, for saving the results. Defaults to a grid of equally spaced values between 0 and max(observed time).

inf_f - inflation factor. If, in order to save computation time, only a sub-sample of the censored observations is used (a case-cohort design), an iflation factor should be provided for the sampled observations (an IPW approach). Defaults to 1 (all censored observations are used).

**The estimation procedure currently does not support ties in the observed event times. If the number of ties in each observed event is small, an approximated solution can be obtained by adding a small amount of noise to break ties (it should be small enough to prevent observed event times from becoming negative).**

## Example for running the algorithm
Suppose that we have a dataset like the pseudo-dataset, simulated using the gamma frailty sampling procedure. In the following lines, the real parameters used in the simulations were given to the function as starting values (for demonstration).

estimate_all(z12=Z12,z13 = Z13,z23 = Z23,obs_times1=V,obs_times2=W,delta1 = delta1,delta2=delta2,
             delta3=delta3,R=R,RecL = 0.05,
             n_Bootstrap = 10,save_address = "C:/Users/User/Desktop",name=paste0("simulation number_",seed),
             initial_pars = c(theta,g12,g13,g23),initial_H = list(H012_real,H013_real,H023_real),check_times = seq(0,0.61,0.01),
             marginal_13_haz = h13_rates, marginal_13_grid = h13_grid_times)

## Sampling Procedures
The codes for sampling data according to gamma/IG/PS frailties are provided. The data are generated assuming piecewise constant hazards with 2 breaking points, as stated in section 5 of our paper.
If the minimal recruitment age is larger than 0, an additonal piece of code is used to generate the "population", and from there extract a grid of marginal hazards for the marginal healthy-death process. 

## UK Biobank Data
The data used for the analysis in the paper are not publicly available. For access to be granted, it is required to submit an application directly to the UKB organization. Additional details regarding the application process, costs, and requirements can be found at the official webpage of the UKbiobank at: 
https://www.ukbiobank.ac.uk/scientists-3/
