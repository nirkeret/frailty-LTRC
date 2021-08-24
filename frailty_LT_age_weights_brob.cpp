#include <stdio.h>
#include <stdlib.h>
/**#include <math.h>**/
#include <Rcpp.h>
#include <cmath>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

//#include "frain.h"

using namespace Rcpp;
using namespace std;
//#define E 2.7182818284590452353602874713527


#define likely(x) __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)

#define Y(times,i,t) times[i] >= t

//#define MY_ADD_LN(x,y) (isinf(x) && isinf(y)? y: (x > y? x + log(1+exp(y-x)): y + log(1 + exp(x-y))))
double MY_ADD_LN(double x,double y) {
	if (isinf(x) && isinf(y)) {
		return y;
	} else {
		if (x > y) {
			return x + log(1+exp(y-x));
		}
		else {
			return y + log(1 + exp(x-y));
		}
	}
}

double* R2C_vec(NumericVector a) {
	int n = a.size();
	double* a_c = (double*)malloc(sizeof(double)*n);
	for(int i = 0; i < n; ++i) {
        	a_c[i] = a[i];
  	}
	return a_c;
}

double* R2C_mat(NumericMatrix a) {
	int n = a.nrow();
	int m = a.ncol();
	double* a_c = (double*)malloc(sizeof(double)*n*m);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
        		a_c[i*m + j] = a(i,j);
		}
  	}
	return a_c;
}


/**#define N(times,delta,i,t) (times[i] < t && delta[i] != 0)? 1:0 
**/

int N(double* times, double* delta, int i, double t) {
	if (times[i] < t and delta[i] > 0) {
		return 1;
	}
	return 0;
}

inline
double exp4(register double x) {
  x = 1.0 + x / 8192;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x;
  return x;
}

#define my_exp(x) exp4(x)


inline int find_m_pos(double f, double* __restrict__ arr, int size) {
	int n = size/2;	
	while (!(arr[n] <= f and f <= arr[n + 1])) {
		if (f > arr[n]) n = size + (size-n) / 2;
		else n = n/2;
	}

	return n;
}


int binarySearch(double* __restrict__ arr, int l, int r, double x)
{
   if (r >= l)
   {
        int mid = l + (r - l)/2;
 
        // If the element is present at the middle itself
        if (arr[mid] == x)  return mid;
 
        // If element is smaller than mid, then it can only be present
        // in left subarray
        if (arr[mid] > x) return binarySearch(arr, l, mid-1, x);
 
        // Else the element can only be present in right subarray
        return binarySearch(arr, mid+1, r, x);
   }
 
   // We reach here when element is not present in array
   return r;
}

/**
int sum_e_matrix(double* __restrict__ a, int n, double* __restrict__  b,int m, double* times, double* __restrict__  sortedTimes, double* __restrict__ res)
{
	sum_e_matrix2(a,n,b,m,times,sortedTimes,res);
	return 0;
}**/





int get_gz_vec(double* __restrict__ z,  double* __restrict__  g, double* __restrict__  res,int n_rows, int n_variables) {
	for (int i=0; i < n_rows; i++) {
		res[i] = 0;
		for (int j=0; j < n_variables; j++) {
			res[i] += z[i*n_variables + j] * g[j];
		}
	}
	return 0;	
}


int get_exp_vec(double* __restrict__ in,  double* __restrict__  out,int len) {
	for (int i=0; i < len; i++) {
		out[i] = exp(in[i]);
	}
	return 0;	
}
/**

int derive_A_combined_c(
			double theta,
			double* __restrict__ z,
			double* __restrict__  g,
			double* __restrict__  H_times,
			double* __restrict__ times,
			double* __restrict__  sortedTimes,
			double* __restrict__ res,
			int n_rows,
			int n_columns,
			int n_variables)
{
	cuda_derive_A_combined_c(theta,
			z,
			g,
			H_times,
			times,
			sortedTimes,
			res,
			n_rows,
			n_columns,
			n_variables);

	return 0;
}
**/


// Derive A matrix by theta and all gamma variales 
// Res expected to be a n*p+1 matrix where p is the number of explenatory variables.
// The first column of res will be the derivative by theta
// The next columns will be the derivative by gamma
int derive_A_combined_c(
			double theta,
			double* __restrict__ z_a,
			double* __restrict__ z_b,
			double* __restrict__  g_a,
			double* __restrict__  g_b,
			double* __restrict__  H_times_a,
			double* __restrict__  H_times_b,
			double* __restrict__  h_times_a,
			double* __restrict__  O_times,
			double* __restrict__  grid_times,
			double* __restrict__ res,
			int n_rows,
			int n_columns,
			int n_variables_a,int n_variables_b)
{
	int n_variables = n_variables_a + n_variables_b;
	for (int i=0; i < n_rows; i++) {
		for (int j=0; j <= n_variables + 1; j++){
			res[i*(n_variables + 2) + j] = -INFINITY;
		}
	}


	double diff[n_columns];
	double gz_vec_a[n_rows];
	double gz_vec_exp_a[n_rows];
	double gz_vec_b[n_rows];
	double gz_vec_exp_b[n_rows];


        for (int i=1; i < n_columns; i++) {
                diff[i] = grid_times[i]-grid_times[i-1];
        }
        diff[0] = grid_times[0];

	get_gz_vec(z_a,g_a,gz_vec_a,n_rows,n_variables_a);
	get_exp_vec(gz_vec_a,gz_vec_exp_a,n_rows);

	get_gz_vec(z_b,g_b,gz_vec_b,n_rows,n_variables_b);
	get_exp_vec(gz_vec_b,gz_vec_exp_b,n_rows);

        #pragma omp parallel for
        for (int i=0; i < n_rows; i++) {
		register double val = 0;
                register double cur_t = O_times[i];
                for (register int j=0;  j< n_columns; j++) {
                        if (unlikely(cur_t <= grid_times[j])) {
				register double chopchick =  j > 0? cur_t-grid_times[j-1]: cur_t;
				val = log(h_times_a[j]*chopchick) + gz_vec_a[i] + theta*(H_times_a[j]*gz_vec_exp_a[i] + H_times_b[j]*gz_vec_exp_b[i]);
				
                        } else {
                        	val = log(h_times_a[j]*diff[j]) + gz_vec_a[i] + theta*(H_times_a[j]*gz_vec_exp_a[i] + H_times_b[j]*gz_vec_exp_b[i]);
				
			}
			
			res[i*(n_variables + 2)] = MY_ADD_LN(res[i*(n_variables + 2)], val);			
			res[i*(n_variables + 2) + 1] = MY_ADD_LN(res[i*(n_variables + 2) + 1], val + log(gz_vec_exp_a[i]*H_times_a[j] + gz_vec_exp_b[i]*H_times_b[j]));
						
			for (register int k=2; k <= n_variables_a + 1; k++) {
				res[i*(n_variables + 2) + k] = MY_ADD_LN(res[i*(n_variables + 2) + k], val + log((abs(z_a[i*n_variables_a + k-2]) + theta*abs(z_a[i*n_variables_a + k-2])*gz_vec_exp_a[i]*H_times_a[j]))); 
			}
			for (register int k= 2 + n_variables_a; k <= n_variables + 1; k++) {
				res[i*(n_variables + 2) + k] = MY_ADD_LN(res[i*(n_variables + 2) + k], val + log((theta*abs(z_b[i*n_variables_b + k-2 - n_variables_a])*gz_vec_exp_b[i]*H_times_b[j]))); 
			}
			
			if (unlikely(cur_t <= grid_times[j])) break;
                }

        }
	return 0;
}

// [[Rcpp::export]]
NumericMatrix derive_A_combined(NumericVector theta_hat, NumericMatrix z_a,NumericMatrix z_b, NumericVector g_a,NumericVector g_b, NumericVector H_times_a,NumericVector H_times_b,NumericVector h_times_a, NumericVector O_times, NumericVector grid_times) {
  int n_rows = z_a.nrow();
  int n_columns = grid_times.size();
  int n_variables_a = z_a.ncol();
  int n_variables_b = z_b.ncol();
  
  double* z_a_c = R2C_mat(z_a);
  double* z_b_c = R2C_mat(z_b);
  double* g_a_c = R2C_vec(g_a);
  double* g_b_c = R2C_vec(g_b);
  double* H_times_a_c = R2C_vec(H_times_a);
  double* H_times_b_c = R2C_vec(H_times_b); 
  double* h_times_a_c =  R2C_vec(h_times_a);
  double* O_times_c =  R2C_vec(O_times);
  double* grid_times_c = R2C_vec(grid_times);
  //double out_c[n_rows][n_variables + 1];
  double out_c[n_rows * (n_variables_a + n_variables_b + 2)];
  double theta = theta_hat[0];

  NumericMatrix out(n_rows,n_variables_a + n_variables_b + 2);

  derive_A_combined_c(
    //derive_A_combined_cuda_int(
	theta,z_a_c,z_b_c,g_a_c,g_b_c,H_times_a_c,H_times_b_c,h_times_a_c,O_times_c,grid_times_c,out_c,n_rows,n_columns,n_variables_a,n_variables_b
	);


  free(z_a_c);
  free(z_b_c);
  free(g_a_c);
  free(g_b_c);
  free(H_times_a_c);
  free(H_times_b_c);
  free(h_times_a_c);
  free(O_times_c);
  free(grid_times_c);
 

  for(int i = 0; i < n_rows; ++i) {
	for(int j = 0; j <= n_variables_a + n_variables_b + 1; ++j) {
    		out(i,j) = out_c[i*(n_variables_a + n_variables_b + 2) + j];
	}
  }
  
  return out;
}




// Derive A matrix by theta and all gamma variales 
// Res expected to be a n*p+1 matrix where p is the number of explenatory variables.
// The first column of res will be the derivative by theta
// The next columns will be the derivative by gamma
int derive_A_combined_23_c(
			double theta,
			double* __restrict__ z,
			double* __restrict__  g,
			double* __restrict__  H_times,
			double* __restrict__  h_times,
			double* __restrict__ O_times,
			double* __restrict__  grid_times,
			double* __restrict__ res,
			int n_rows,
			int n_columns,
			int n_variables)
{
	for (int i=0; i < n_rows; i++) {
		for (int j=0; j <= n_variables + 2; j++){
			res[i*(n_variables + 3) + j] = -INFINITY;
		}
	}


	double diff[n_columns];
	double gz_vec[n_rows];
	double gz_vec_exp[n_rows];

        for (int i=1; i < n_columns; i++) {
                diff[i] = grid_times[i]-grid_times[i-1];
        }
        diff[0] = grid_times[0];

	get_gz_vec(z,g,gz_vec,n_rows,n_variables);
	get_exp_vec(gz_vec,gz_vec_exp,n_rows);


        #pragma omp parallel for
        for (int i=0; i < n_rows; i++) {
		register double val = 0;
                register double cur_t = O_times[i];
                for (register int j=0;  j< n_columns; j++) {
                        if (unlikely(cur_t <= grid_times[j])) {
				register double chopchick =  j > 0? cur_t-grid_times[j-1]: cur_t;
				val = log(h_times[j]*chopchick) + gz_vec[i] + (theta/(1 + theta))*gz_vec_exp[i]*H_times[j] - log(1 + theta);
                        } else {
                        	val = log(h_times[j]*diff[j]) + gz_vec[i] + (theta/(1 + theta))*gz_vec_exp[i]*H_times[j] - log(1 + theta);
			}
			res[i*(n_variables + 3)] = MY_ADD_LN(res[i*(n_variables + 3)], val);
			//res[i*(n_variables + 3) + 1] = MY_ADD_LN(res[i*(n_variables + 3) + 1], val + -2*log(1+);
			if ((gz_vec_exp[i]*H_times[j]/(1+theta) - 1) > 0) {
				res[i*(n_variables + 3) + 1] = MY_ADD_LN(res[i*(n_variables + 3) + 1], val -log(1+theta) + log(gz_vec_exp[i]*H_times[j]/(1+theta) - 1));
			} else {
				res[i*(n_variables + 3) + 2] = MY_ADD_LN(res[i*(n_variables + 3) + 2], val -log(1+theta) + log(1- (gz_vec_exp[i]*H_times[j]/(1+theta))));
			}
			for (register int k=3; k <= n_variables + 2; k++) {
				res[i*(n_variables + 3) + k] = MY_ADD_LN(res[i*(n_variables + 3) + k], val + log((abs(z[i*n_variables + k-3]) + (theta/(1 + theta))*abs(z[i*n_variables + k-3])*gz_vec_exp[i]*H_times[j]))); 
			}
			if (unlikely(cur_t <= grid_times[j])) break;
                }
        }

	return 0;
}

// [[Rcpp::export]]
NumericMatrix derive_A_combined_23(NumericVector theta_hat, NumericMatrix z, NumericVector g, NumericVector H_times, NumericVector h_times, NumericVector O_times, NumericVector grid_times) {
  int n_rows = z.nrow();
  int n_columns = grid_times.size();
  int n_variables = z.ncol();
  int tmp = g.size();
  //printf("*** %d\n", tmp);

  double* z_c = R2C_mat(z);
  double* g_c = R2C_vec(g);
 
  double* H_times_c = R2C_vec(H_times);
  double* h_times_c = R2C_vec(h_times);
  double* O_times_c =  R2C_vec(O_times);
  double* grid_times_c = R2C_vec(grid_times);
  //double out_c[n_rows][n_variables + 1];
  double out_c[n_rows * (n_variables + 3)];
  double theta = theta_hat[0];

  NumericMatrix out(n_rows,n_variables + 3);
  //printf("1\n");
    derive_A_combined_23_c(
    //derive_A_combined_23_cuda_int(
		theta,z_c,g_c,H_times_c,h_times_c,O_times_c,grid_times_c,(double*)out_c,n_rows,n_columns,n_variables);
  //printf("2\n");

  free(z_c);
  free(g_c);
  free(H_times_c);
  free(h_times_c);
  free(O_times_c);
  free(grid_times_c);

  for(int i = 0; i < n_rows; ++i) {
	for(int j = 0; j <= n_variables + 2; ++j) {
    		out(i,j) = out_c[i*(n_variables + 3) + j];
	}
  }

  return out;
}




inline void diff_vec(double* __restrict__ vec , double* __restrict__ ret, int size) {
	for (int i = 0; i < size; i++) {
		ret[i] = vec[i + 1] - vec[i];
	}
}

inline void diff_vec_compete(double* __restrict__ vec , double* __restrict__ ret, int size) {
	for (int i = 0; i < size - 1; i++) {
		ret[i] = vec[i+1] - vec[i];
	}
}

inline void copy_num_vector(NumericVector vec , double* __restrict__ ret) {
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		ret[i] = vec[i];
	}
}



#define MAX(a,b) ((a) > (b) ? a : b)







void get_ln_A_col_c(double* __restrict__ egz_a,
		  double* __restrict__ egz_b,
		  double theta_hat,
		  double* __restrict__ O_time,
		  double* __restrict__ h_a,
		  double* __restrict__ H_a,
		  double* __restrict__ H_b,
		  double* __restrict__ grid_times,
		  int N, int N_E,double* res) {
	for (int i = 0; i < N; i++) {
		res[i] = 0;
	}
	
	double diff[N_E];

        for (int i=1; i < N_E; i++) {
                diff[i] = grid_times[i]-grid_times[i-1];
        }
        diff[0] = grid_times[0];

        for (int sub_indx=0; sub_indx < N; sub_indx++) {
                register double sum = 0;
                register double cur_t = O_time[sub_indx];
                if (unlikely(cur_t <= grid_times[0])) {
                        sum += cur_t;
                } else {
                        sum += diff[0];
                }
		sum *= h_a[0]*egz_a[sub_indx]*exp(theta_hat*(egz_a[sub_indx]*H_a[0] + egz_b[sub_indx]*H_b[0]));
                res[sub_indx] = log(sum);
        }

        #pragma omp parallel for
        for (int sub_indx=0; sub_indx < N; sub_indx++) {
                register double sum = res[sub_indx];
                register double cur_t = O_time[sub_indx];
                if (cur_t <= grid_times[0]) {
                        continue;
                }
                for (register int event_indx=1;  event_indx< N_E; event_indx++) {
                      if (unlikely(cur_t < grid_times[event_indx])) {
				register double tmp = log(
						    (cur_t - grid_times[event_indx-1])*h_a[event_indx]*egz_a[sub_indx]) +
						    theta_hat*(
								egz_a[sub_indx]*H_a[event_indx] + 
								egz_b[sub_indx]*H_b[event_indx]);
                         	if (sum > tmp) {
			 		sum = sum + log(1 + exp(tmp-sum));
			 	}
				else {
					sum = tmp + log(1 + exp(sum-tmp));
			 	}
				break;
                        }
			register double tmp = log(diff[event_indx]*h_a[event_indx]*egz_a[sub_indx]) + theta_hat*(egz_a[sub_indx]*H_a[event_indx] + egz_b[sub_indx]*H_b[event_indx]);
			if (sum > tmp) {
				sum = sum + log(1 + exp(tmp-sum));
			}
			else {
				sum = tmp + log(1 + exp(sum-tmp));
			}
                }
                res[sub_indx] = sum;
        }
}



// [[Rcpp::export]]
NumericVector get_ln_A_col(NumericVector egz_a, NumericVector egz_b, NumericVector theta, NumericVector O_times,  NumericVector h_a,  NumericVector H_a,NumericVector H_b ,NumericVector gridtimes) {
  int N = egz_a.size();
  int N_E = gridtimes.size();

  double* egz_a_c = R2C_vec(egz_a);
  double* egz_b_c = R2C_vec(egz_b);
  double thetha =  theta[0];
  double* O_times_c = R2C_vec(O_times);
  double* h_a_c = R2C_vec(h_a);
  double* H_a_c = R2C_vec(H_a);
  double* H_b_c = R2C_vec(H_b);
  double* gridtimes_c = R2C_vec(gridtimes);
  double* out_c = (double* )malloc(sizeof(double)*N);

  NumericVector out(N);

  get_ln_A_col_c( 
   //get_ln_A_col_cuda_int(
		egz_a_c,
		egz_b_c,
	        thetha,
		O_times_c,
		h_a_c,
		H_a_c,
		H_b_c,
		gridtimes_c,
		N,N_E,out_c);



  free(egz_a_c);
  free(egz_b_c);
  free(O_times_c);
  free(h_a_c);
  free(H_b_c);
  free(gridtimes_c);


  for(int i = 0; i < N; ++i) {
    out[i] = out_c[i];
  }

  free(out_c);

  return out;
}


void get_ln_A_col_23(double* __restrict__ egz,
		  double theta_hat,
		  double* __restrict__ O_time,
		  double* __restrict__ h,
		  double* __restrict__ H,
		  double* __restrict__ grid_times,
		  int N, int N_E,double* res) {
	for (int i = 0; i < N; i++) {
		res[i] = 0;
	}
	
	double diff[N_E];

        for (int i=1; i < N_E; i++) {
                diff[i] = grid_times[i]-grid_times[i-1];
        }
        diff[0] = grid_times[0];

        for (int sub_indx=0; sub_indx < N; sub_indx++) {
                register double sum = 0;
                register double cur_t = O_time[sub_indx];
                if (unlikely(cur_t < grid_times[0])) {
                        sum += cur_t;
                } else {
                        sum += diff[0];
                }
		sum *= (1/(1+theta_hat))*h[0]*egz[sub_indx]*exp(egz[sub_indx]*H[0]*(theta_hat/(1+theta_hat)));
                res[sub_indx] = log(sum);
        }

        #pragma omp parallel for
        for (int sub_indx=0; sub_indx < N; sub_indx++) {
                register double sum = res[sub_indx];
                register double cur_t = O_time[sub_indx];
                if (cur_t <= grid_times[0]) {
                        continue;
                }
                for (register int event_indx=1;  event_indx< N_E; event_indx++) {
                      if (unlikely(cur_t < grid_times[event_indx])) {
				register double tmp = log((cur_t - grid_times[event_indx-1])*(1/(1+theta_hat))*h[event_indx]*egz[sub_indx]) +
							 (egz[sub_indx]*H[event_indx])*(theta_hat/(1+theta_hat));
                         	if (sum > tmp) {
			 		sum = sum + log(1 + exp(tmp-sum));
			 	}
				else {
					sum = tmp + log(1 + exp(sum-tmp));
			 	}
			 	break;
                        }
			register double tmp = log(diff[event_indx]*(1/(1+theta_hat))*h[event_indx]*egz[sub_indx]) +
							 (egz[sub_indx]*H[event_indx])*(theta_hat/(1+theta_hat));
			if (sum > tmp) {
				sum = sum + log(1 + exp(tmp-sum));
			}
			else {
				sum = tmp + log(1 + exp(sum-tmp));
			}

                }
                res[sub_indx] = sum;
        }
}


// [[Rcpp::export]]
NumericVector get_ln_A_col_23(NumericVector egz, NumericVector theta, NumericVector O_times,  NumericVector h,  NumericVector H ,NumericVector gridtimes) {
  int N = egz.size();
  int N_E = gridtimes.size();

  double* egz_c = R2C_vec(egz);
  double thetha =  theta[0];
  double* O_times_c = R2C_vec(O_times);
  double* h_c = R2C_vec(h);
  double* H_c = R2C_vec(H);
  double* gridtimes_c = R2C_vec(gridtimes);
  double* out_c = (double* )malloc(sizeof(double)*N);

  NumericVector out(N);

  get_ln_A_col_23(
    //get_ln_A_col_23_cuda_int(
		egz_c,
	        thetha,
		O_times_c,
		h_c,
		H_c,
		gridtimes_c,
		N,N_E,out_c);



  free(egz_c);
  free(O_times_c);
  free(h_c);
  free(H_c);
  free(gridtimes_c);


  for(int i = 0; i < N; ++i) {
    out[i] = out_c[i];
  }

  free(out_c);

  return out;
}






void c_estimate_H_JK(     double* __restrict__ gz_12,
			  double* __restrict__ gz_13,
                          double theta_hat,
                          double* __restrict__ H_12,
			  double* __restrict__ H_13,
			  double* __restrict__ h_12,
			  double* __restrict__ h_13,
			  double* __restrict__ is_12,
			  double* __restrict__ O_times,
			  double* __restrict__ delta,
			  double* __restrict__ grid_times,
			  double* trunc_subjects,
			  double* weights_subjects,
			  double* weights_events_12,
			  double* weights_events_13,
			  int N_subj,
			  int N_EV,
			  int N_EV_12,
			  int N_EV_13,
			  double* __restrict__ ret_12,double* __restrict__ ret_13) {

	double A_12[N_subj];
	double A_13[N_subj];
	double gz_vec_exp_12[N_subj];
	double gz_vec_exp_13[N_subj];

	double diff[N_EV];

	for (int i = 0; i < N_subj; i++) {
		A_12[i] = 0;
		A_13[i] = 0;
	}

        for (int i=1; i < N_EV; i++) {
                diff[i] = grid_times[i]-grid_times[i-1];
        }
        diff[0] = grid_times[0];


	double E_val = 0;
	register int event_itr = 0;
	register int subjet_itr = 0;

	get_exp_vec(gz_12,gz_vec_exp_12,N_subj);
	get_exp_vec(gz_13,gz_vec_exp_13,N_subj);


	for (event_itr = 0; event_itr < N_EV_12; event_itr++ ) {
		ret_12[event_itr] = 0;

	}

	for (event_itr = 0; event_itr < N_EV_13; event_itr++ ) {
		ret_13[event_itr] = 0;
	}

	//In order to work in parallel we calculate each row of A and add the contributon of this sample to ret (we sum)
	//But we must keep diffrent vectors for ret for diffrent threads so as to not cause collison vetween the threads when they try to update the ret of the same even time
	//So we create a new array(only if using parralization) the size of events*number of threads each thread will work on a shifted temp array
	//In the end a single thread will merge them.
	#if defined(_OPENMP)
	   double* S_private_12 = NULL;
	   double* S_private_13 = NULL;	
	#endif	
	
	#pragma omp parallel private(event_itr,subjet_itr,E_val)
	{
		#if defined(_OPENMP)
		    const int nthreads = omp_get_num_threads();
		    const int ithread = omp_get_thread_num();

		    #pragma omp single 
		    {
			S_private_12 = new double[N_EV_12*nthreads];
			S_private_13 = new double[N_EV_13*nthreads];
			for(int i=0; i<(N_EV_12*nthreads); i++) S_private_12[i] = 0;
			for(int i=0; i<(N_EV_13*nthreads); i++) S_private_13[i] = 0;
		    }
		#endif	


		#pragma omp for 
		for (subjet_itr = 0; subjet_itr < N_subj; subjet_itr++) {
			int e12_count = 0;
			for (event_itr = 0; event_itr < N_EV; event_itr++ ) {
				register double tmp = exp(theta_hat*(gz_vec_exp_12[subjet_itr]*H_12[event_itr] + gz_vec_exp_13[subjet_itr]*H_13[event_itr]));
				register double alpha_12 = h_12[event_itr]*tmp*gz_vec_exp_12[subjet_itr];
				register double alpha_13 = h_13[event_itr]*tmp*gz_vec_exp_13[subjet_itr];
				A_12[subjet_itr] += alpha_12*diff[event_itr];
				A_13[subjet_itr] += alpha_13*diff[event_itr];
				if (Y(O_times,subjet_itr,grid_times[event_itr])) {
					if (grid_times[event_itr] >= trunc_subjects[subjet_itr]) {
						E_val = (1.0/theta_hat + N(O_times, delta ,subjet_itr, grid_times[event_itr])) / 
						        (1.0/theta_hat + A_12[subjet_itr] + A_13[subjet_itr]);
						if (is_12[event_itr]) {
							double astar = gz_vec_exp_12[subjet_itr] * tmp;
							//if (e12_count == 5) {
							//	printf("%d %d: %f %f %f %f %f\n",event_itr, subjet_itr, E_val,astar,A_12[subjet_itr] , A_13[subjet_itr],weights_subjects[subjet_itr]);
							//}
							#if defined(_OPENMP)		
							S_private_12[ithread*N_EV_12 + e12_count] += E_val*astar*weights_subjects[subjet_itr];
							#else
							//printf("%d\n",e12_count);
							ret_12[e12_count] += E_val*astar*weights_subjects[subjet_itr];
							#endif
		
						} else {
							double astar = gz_vec_exp_13[subjet_itr] * tmp;

							#if defined(_OPENMP)		
							S_private_13[ithread*N_EV_13 + event_itr - e12_count] += E_val*astar*weights_subjects[subjet_itr];
							#else
							ret_13[event_itr - e12_count] += E_val*astar*weights_subjects[subjet_itr];
							#endif							
						}
					}			
				}
				else {
					break;
				}
				if (is_12[event_itr]) {
					e12_count += 1;	
					//if (subjet_itr == 80) {
					//	printf("%d %d:%f %f\n",event_itr, subjet_itr, H_13[event_itr], H_12[event_itr]);
					//}
				}			
			}
		}

	    #if defined(_OPENMP)
	    #pragma omp single
	    {	
	    	for(int i=0; i<N_EV_12; i++) {
			for(int t=0; t<nthreads; t++) {
			    ret_12[i] += S_private_12[N_EV_12*t + i];
			}
	    	}
	    	delete S_private_12;
	    	for(int i=0; i<N_EV_13; i++) {
			for(int t=0; t<nthreads; t++) {
			    ret_13[i] += S_private_13[N_EV_13*t + i];
			}
	    	}
	    	delete S_private_13;
	    }
	    #endif

	}

	for (event_itr = 0; event_itr < N_EV_12; event_itr++ ) {
		ret_12[event_itr] = weights_events_12[event_itr] / ret_12[event_itr];
	}
	for (event_itr = 0; event_itr < N_EV_13; event_itr++ ) {
		ret_13[event_itr] = weights_events_13[event_itr] / ret_13[event_itr];
	}
}


// [[Rcpp::export]]
Rcpp::List estimate_H_JK(NumericVector gz_12,
			      NumericVector gz_13,
			      NumericVector theta, 
			      NumericVector H_12,
			      NumericVector H_13,
			      NumericVector h_12,
			      NumericVector h_13,
			      NumericVector is_12,
			      NumericVector O_times,
			      NumericVector delta,
			      NumericVector grid_times,
			      NumericVector trunc_subjects,
			      NumericVector weights_subjects,
			      NumericVector weights_events_12,
			      NumericVector weights_events_13
) {
  int N = gz_12.size();
  int N_E = grid_times.size();

  int N_E_12 = 0;
  for (int i = 0; i < N_E; i++) {
	if (is_12[i]) {
		N_E_12 += 1;
	}
  }
  int N_E_13 = N_E - N_E_12;

  double* gz_12_c = R2C_vec(gz_12);
  double* gz_13_c = R2C_vec(gz_13);
  double theta_hat =  theta[0];
  double* H_12_c = R2C_vec(H_12);
  double* H_13_c = R2C_vec(H_13);
  double* h_12_c = R2C_vec(h_12);
  double* h_13_c = R2C_vec(h_13);
  double* is_12_c = R2C_vec(is_12);
  double* O_times_c = R2C_vec(O_times);
  double* delta_c = R2C_vec(delta);
  double* grid_times_c = R2C_vec(grid_times);
  double* trunc_subjects_c = R2C_vec(trunc_subjects);
  double* weights_subjects_c = R2C_vec(weights_subjects);
  double* weights_events_12_c = R2C_vec(weights_events_12);
  double* weights_events_13_c = R2C_vec(weights_events_13);

  double* out_12_c = (double* )malloc(sizeof(double)*N_E_12);
  double* out_13_c = (double* )malloc(sizeof(double)*N_E_13);

  NumericVector out12(N_E_12);
  NumericVector out13(N_E_13);

  c_estimate_H_JK( gz_12_c,gz_13_c,
	        theta_hat,
		H_12_c,H_13_c,
		h_12_c,h_13_c,is_12_c,O_times_c,delta_c,
		grid_times_c,trunc_subjects_c,weights_subjects_c,weights_events_12_c,weights_events_13_c,
		N,N_E,N_E_12,N_E_13,out_12_c,out_13_c);


  free(gz_12_c);
  free(gz_13_c);
  free(H_12_c);
  free(H_13_c);
  free(h_12_c);
  free(h_13_c);
  free(is_12_c);
  free(O_times_c);
  free(delta_c);
  free(grid_times_c);
  free(trunc_subjects_c);
  free(weights_subjects_c);
  free(weights_events_12_c);
  free(weights_events_13_c);


  for(int i = 0; i < N_E_12; ++i) {
    out12[i] = out_12_c[i];
  }

  for(int i = 0; i < N_E_13; ++i) {
    out13[i] = out_13_c[i];
  }
  free(out_12_c);
  free(out_13_c);
  return Rcpp::List::create(Rcpp::Named("h12") = out12,
                          Rcpp::Named("h13") = out13);
}





void c_estimate_H_23(double* __restrict__ gz,
                          double theta_hat,
                          double* __restrict__ H,
			  double* __restrict__ grid_times,
			  double* __restrict__ O_times,
			  double* __restrict__ V_times,
			  double* __restrict__ delta,
			  double* __restrict__ A_1_tau,
			  double* trunc_subjects,
			  double* weights_subjects,
			  double* weights_events,
			  int N_subj,
			  int N_EV,
			  double* __restrict__ ret) {

	double a_star[N_subj];
	double A[N_subj];
	double gz_vec_exp[N_subj];
	double h[N_EV];
		
	double E_val = 0;
	int event_itr = 0;
	int subjet_itr = 0;

	diff_vec(H,h,N_EV);
	get_exp_vec(gz,gz_vec_exp,N_subj);
	
	for (event_itr = 0; event_itr < N_EV; event_itr++ ) {
		ret[event_itr] = 0;
	}


	ret[0] = 0;
	for (subjet_itr = 0; subjet_itr < N_subj; subjet_itr++) {
		A[subjet_itr] = A_1_tau[subjet_itr];
	}
	
	//In order to work in parallel we calculate each row of A and add the contributon of this sample to ret (we sum)
	//But we must keep diffrent vectors for ret for diffrent threads so as to not cause collison vetween the threads when they try to update the ret of the same even time
	//So we create a new array(only if using parralization) the size of events*number of threads each thread will work on a shifted temp array
	//In the end a single thread will merge them.
	#if defined(_OPENMP)
	   double* S_private = NULL;	
	#endif	
	
	#pragma omp parallel private(event_itr,subjet_itr,E_val)
	{
		#if defined(_OPENMP)
		    const int nthreads = omp_get_num_threads();
		    const int ithread = omp_get_thread_num();

		    #pragma omp single 
		    {
			S_private = new double[N_EV*nthreads];
			for(int i=0; i<(N_EV*nthreads); i++) S_private[i] = 0;
		    }
		#endif	

		#pragma omp for 
		for (subjet_itr = 0; subjet_itr < N_subj; subjet_itr++) {
			for (event_itr = 0; event_itr < N_EV; event_itr++ ) {
				register double astar = (1/(1+ theta_hat))*gz_vec_exp[subjet_itr]*exp(gz_vec_exp[subjet_itr]*H[event_itr]*(theta_hat/(1+ theta_hat)));
				if (Y(O_times,subjet_itr,grid_times[event_itr])) {	
					A[subjet_itr] += astar * h[event_itr];
				}
				else {
					break;
				}
				if (Y(O_times,subjet_itr,grid_times[event_itr]) && 
				   (grid_times[event_itr] > V_times[subjet_itr]) &&
				   (V_times[subjet_itr] >= trunc_subjects[subjet_itr] || ((trunc_subjects[subjet_itr]) < grid_times[event_itr]))) {
					E_val = (1.0/theta_hat + 1 + N(O_times, delta ,subjet_itr, grid_times[event_itr])) / (1.0/theta_hat + A[subjet_itr]);
					#if defined(_OPENMP)		
					S_private[ithread*N_EV +event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#else
					ret[event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#endif	
				}
			}
		}


	
	    #if defined(_OPENMP)
	    #pragma omp single
	    {	
	    	for(int i=0; i<N_EV; i++) {
			for(int t=0; t<nthreads; t++) {
			    ret[i] += S_private[N_EV*t + i];
			}
	    	}
	    	delete S_private;
	    }
	    #endif

	}

	for (event_itr = 0; event_itr < N_EV; event_itr++ ) {
		ret[event_itr] = weights_events[event_itr] / ret[event_itr];
	}
}



// [[Rcpp::export]]
NumericVector estimate_H_23(NumericVector gz,
			      NumericVector theta, 
			      NumericVector H,
			      NumericVector grid_times,
			      NumericVector O_times,
			      NumericVector V_times,
			      NumericVector delta,
			      NumericVector A_1_tau,
			      NumericVector trunc_subjects,
			      NumericVector weights_subjects,
			      NumericVector weights_events
) {
  int N = gz.size();
  int N_E = grid_times.size();


  double* gz_c = R2C_vec(gz);
  double theta_hat =  theta[0];
  double* H_c = R2C_vec(H);
  double* grid_times_c = R2C_vec(grid_times);
  double* O_times_c = R2C_vec(O_times);
  double* V_times_c = R2C_vec(V_times);
  double* delta_c = R2C_vec(delta);
  double* A_1_tau_c = R2C_vec(A_1_tau);
  double* trunc_subjects_c = R2C_vec(trunc_subjects);
  double* weights_subjects_c = R2C_vec(weights_subjects);
  double* weights_events_c = R2C_vec(weights_events);

	
  double* out_c = (double* )malloc(sizeof(double)*N_E);

  NumericVector out(N_E);

  c_estimate_H_23(gz_c,
	        theta_hat,
		H_c,
		grid_times_c,
		O_times_c,V_times_c,delta_c,A_1_tau_c,trunc_subjects_c,weights_subjects_c,weights_events_c,
		N,N_E,out_c);


  free(gz_c);
  free(H_c);
  free(grid_times_c);
  free(O_times_c);
  free(V_times_c);
  free(delta_c);
  free(A_1_tau_c);
  free(trunc_subjects_c);
  free(weights_subjects_c);
  free(weights_events_c);



  for(int i = 0; i < N_E; ++i) {
    out[i] = out_c[i];
  }

  free(out_c);
  return out;
}



