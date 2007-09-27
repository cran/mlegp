/***************************************************************************
    Maximum Likelihood Estimation of Gaussian Processes (mlegp)
    Copyright (C) 2007, Garrett Dancik

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    For questions or comments, contact Garrett Dancik (gdancik@iastate.edu)
****************************************************************************/

/***************************************************************************
fit_gp.h - fits a gp to a set of observations, using numerical methods and 
	closed form solutions to obtain mles
***************************************************************************/

#ifndef __FIT_GP__
#define __FIT_GP__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

#include <time.h>
#include <math.h>

#ifdef __useR__
#include <R_ext/Utils.h>
#endif

#include "matrix_vector.h"
#include "gp.h"
#include "print.h"


// parameters that do not vary; will be passed into minimization functions (f_min, etc)
typedef struct {
	gsl_matrix *X;       	// the design that produces the observations
	gsl_matrix *fX;      	// the mean function will have the form (fX)'B
	gsl_matrix *Y;		// the observations
	int num_params;		// number of columns of fX 
	double h;		// used in derivative approximation [f(x+h) - f(x+h)] / h
	double *nugget;	        // pointer to nugget 'matrix'; NULL otherwise
	int verbose;		// level of detail for output
	int constantNugget;     // 1 if we are estimating a constant nugget	
	double min_nugget;	// minimum value of the nugget, or can be used to force a nugget
} gp_params; 


/* vector_exp_check - takes exp of each element in the vector; sets to 0 if element < -500 */
void vector_exp_check(gsl_vector *v) {
	int i = 0;
	for (i = 0; i < v->size; i++) {
		if (v->data[i] < -500) v->data[i] = 0;
		else v->data[i] = exp(v->data[i]);
	}
}

/* vector_log_subset - takes log of elements 0 through elements (size-1) in vector v */
void vector_log_subset(gsl_vector *v, int size) {
	int i = 0; 
	for (i = 0; i < size; i++) {
		v->data[i] = log(v->data[i]);
	}
}



/**********************************************************************************
f_min - finds the neg log like of the gp
	gsl_vector *orig_v - vector of parameters being estimated, consisting of
		correlation coefficients, overall sig2, sig2(nugget), in that
		order; all parameters must be on log scale; sig2(nuget) is 
		a scaled nugget and is actually sig2(nugget) / sig2(gp) 
	void *params - pointer to gp_params structure, which determines whether
			we are using a constant nugget, nugget matrix, the
			mean function, etc.
	
	returns: -log likelihood of gp; or DBL_MAX if there is an error inverting 
			var-cov matrix
***********************************************************************************/
	
double f_min(const gsl_vector *orig_v, void *params) {
	int i = 0;
	gp_params *p = (gp_params *) params;

	gsl_matrix *corr = gsl_matrix_alloc(p->X->size1,p->X->size1);

	gsl_vector *v = gsl_vector_alloc(orig_v->size);
	for (i = 0; i < v->size; i++) {
		v->data[i] = orig_v->data[i];
	}
	
	vector_exp_check(v);   // transform back to original values
	gsl_vector *B = gsl_vector_alloc(p->num_params);
	for (i = 0; i < B->size; i++) {
		B->data[i] = v->data[i];
	}

	createCorrMatrix(p->X, B, corr);

	if (p->nugget != NULL) addNuggetMatrix(corr, v->data[p->num_params], p->nugget);
	else if (p->constantNugget != 1) addNugget(corr, v->data[p->num_params]);
	addNugget(corr, p->min_nugget); // add the minimum nugget

	gsl_matrix *corrInv = gsl_matrix_alloc(p->X->size1,p->X->size1);
	
	gsl_set_error_handler_off();

	i = solve(corr, corrInv);
	
	gsl_set_error_handler(NULL);

	if (i == 1) {  // singular matrix!
		gsl_vector_free(B);
		gsl_matrix_free(corr);
		gsl_matrix_free(corrInv);
		return DBL_MAX;
	}
		
	gsl_matrix *bhat = gsl_matrix_alloc(p->fX->size2, 1);
	
	if (calcBhat(p->fX, corrInv, p->Y, bhat) == 1) {
		gsl_vector_free(B);
		gsl_matrix_free(corr);
		gsl_matrix_free(corrInv);
		return DBL_MAX;
	}

	gsl_matrix *mu_hat = gsl_matrix_alloc(p->fX->size1, bhat->size2);
	matrix_multiply(p->fX, bhat, mu_hat);


	double sig2 = calcMLESig2(p->Y, mu_hat, corrInv);

	gsl_matrix_scale(corr, sig2);

	gsl_vector *Y_vec = gsl_vector_alloc(p->Y->size1);
	gsl_vector *mu_hat_vec = gsl_vector_alloc(p->Y->size1);

	toVector(p->Y, Y_vec);
	toVector(mu_hat, mu_hat_vec);

	double ans = - logdmvnorm(Y_vec, mu_hat_vec, corr);
	gsl_matrix_free(corr);
	gsl_vector_free(B);
	gsl_vector_free(v);
	gsl_matrix_free(corrInv);
	gsl_matrix_free(bhat);
	gsl_matrix_free(mu_hat);
	gsl_vector_free(Y_vec);
	gsl_vector_free(mu_hat_vec);
	return ans;

} // end f_min




// finite difference approximation to gradient; uses forward difference
// if current v has singular var-cov matrix, set gradient to 0;
// if f(v+h) has singular var-cov matrix, try backward difference approx,
// if f(v-h) has singular var-cov matrix, set gradient to 0

/***********************************************************************************
df_min - approximates partial derivatives using a forward difference approximation
	for each parameter using the current values of parameters given in the 
	vecotr 'v'. These gradients are returned in the vector df; if forward
	approximation does not work; try backward difference approximation. If 
	none of these work or current v cannot be evaluated, return a gradient 
	0. This function is not called until after simplex method, and setting
	the gradient to 0 is consistent with having an mle near the edge of 
	a parameter space where the var-cov matrix becomes singular.

	gsl_vector *v and void *params are the same as in f_min
************************************************************************************/

void df_min(const gsl_vector *v, void *params, gsl_vector *df) {
	gsl_vector *v_plus_h = gsl_vector_alloc(v->size);
	gp_params *p = (gp_params *) params;

	int i;
	double fv_plus_h, fv;
	for (i = 0; i < v->size; i++) {
		gsl_vector_memcpy(v_plus_h, v);
		//v_plus_h->data[i] = v->data[i] + p->h;        // when v is on actual scale
		// v contains log of all parameters
		v_plus_h->data[i] = v->data[i] + log(1 + p->h*exp(-v->data[i])); //we will look at v, v+h on true scale
		fv_plus_h = f_min(v_plus_h, params);
		fv = f_min(v, params); 
		if (fv == DBL_MAX) {
			df->data[i] = 0;
		}
		else if (fv_plus_h == DBL_MAX) {     // try backward approximation
			gsl_vector_memcpy(v_plus_h, v);
			//we will look at v, v-h on true scale
			v_plus_h->data[i] = v->data[i] - log(1 + p->h*exp(-v->data[i])); 
			fv_plus_h = f_min(v_plus_h, params);
			if (fv_plus_h == DBL_MAX) {
			    df->data[i] = 0;
			}
			else df->data[i] = (fv - fv_plus_h) / -p->h;		

		}
		else df->data[i] = (fv_plus_h - fv) / p->h;
	}
}


void fdf_min(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
	*f = f_min(v, params);
	df_min(v, params, df);
}


/****************************************************************************************************************
optimizeUsingSimplex - uses Nelder Mead Simplex to find mles
	gsl_vector *v - vector of parameters, same as in f_min
	gp_params *params - same as in f_min
	gsl_vector *steps - pointer to vector of initial step sizes
	int max_iter - maximum number of iterations
	double break_size - stoppint criteria

	returns: the value at the minimum (neg log likelihood)
****************************************************************************************************************/

double optimizeUsingSimplex(gsl_vector *v, gp_params *params, gsl_vector *steps, int max_iter, double break_size) {

	int every = 20; // display status every 'every' iterations

	/*** first use a simplex approach ***/

	gsl_multimin_function f_minimize;  // simplex
	f_minimize.f = f_min;
	f_minimize.n = v->size;
	f_minimize.params = (void *) params;

	const gsl_multimin_fminimizer_type *min_type = gsl_multimin_fminimizer_nmsimplex;   // simplex

	gsl_multimin_fminimizer *the_minimizer;
	the_minimizer = gsl_multimin_fminimizer_alloc(min_type,v->size);

	gsl_multimin_fminimizer_set(the_minimizer, &f_minimize, v,steps); // simplex
	int iter = 0, status;
	do {
		#ifdef __useR__
		R_CheckUserInterrupt();
		#endif
		iter++;
		status = gsl_multimin_fminimizer_iterate(the_minimizer); // simplex
		if (params->verbose > 0 && iter % every == 0) {
			printout("\titeration: %d, ", iter);
			printout("loglike = %f", -the_minimizer->fval); // simplex 
			if (params->verbose == 2) {
				printout(", est: "); printVector(the_minimizer->x); 
			}
			else printout("\n");
		}
		
		if (status) {
			break;
		}

		status = gsl_multimin_test_size(the_minimizer->size, break_size); // simplex
		
		if (status == GSL_SUCCESS) {
			if (params->verbose > 0) {
				printout("\tminimum found at iter %d, loglike = %f\n", iter, -the_minimizer->fval);
				//printVector(the_minimizer->x);
			}
		}
	
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_vector_memcpy(v, the_minimizer->x);
	double fval = the_minimizer->fval;
	gsl_multimin_fminimizer_free(the_minimizer); // simplex
	return fval;
}

/****************************************************************************************************************
optimizeUsingGradient - uses bfgs method to find mles
	gsl_vector *v - vector of parameters, same as in f_min
	gp_params *params - same as in f_min
	double step_size - size of the first trial step 
	double tol - stopping condition is when dot product of gradient < tol
	int max_iter - maximum number of iterations

	returns: the value at the minimum (neg log likelihood)
****************************************************************************************************************/
double optimizeUsingGradient(gsl_vector *v, gp_params *params, double step_size, double tol, int max_iter) {

	int every = 5;     // report status every so many iterations

	gsl_multimin_function_fdf f_minimize;  
	f_minimize.f = f_min;
	f_minimize.n = v->size;

	f_minimize.params = (void *) params;
	f_minimize.df = df_min;		// gradient only
	f_minimize.fdf = fdf_min;	// gradient only

	const gsl_multimin_fdfminimizer_type *min_type;
	min_type = gsl_multimin_fdfminimizer_vector_bfgs; 
	///min_type = gsl_multimin_fdfminimizer_conjugate_pr; 

	gsl_multimin_fdfminimizer *the_minimizer;
	the_minimizer = gsl_multimin_fdfminimizer_alloc(min_type,v->size);

	gsl_multimin_fdfminimizer_set(the_minimizer, &f_minimize, v,step_size, tol); // gradient

	int iter = 0, status = 0;

	if (max_iter > 0) {
	do {
		#ifdef __useR__
		R_CheckUserInterrupt();
		#endif
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(the_minimizer); // gradient
		if (params->verbose > 0 && (iter % every) == 0) { 
			printout("\titer: %d, loglike = %f", iter, -the_minimizer->f); // gradient
			if (params->verbose == 2) {
				printout(", est: ");
				printVector(the_minimizer->x);
			}
			else printout("\n");
		}
		if (status) {
			break;
		}

		status = gsl_multimin_test_gradient(the_minimizer->gradient, tol); // gradient
		
		if (status == GSL_SUCCESS) {
			if (params->verbose > 0) {
				printout("BFGS method complete, loglike = ");
				printout("%f\n", -the_minimizer->f);
			}
		}
	
	} while (status == GSL_CONTINUE && iter < max_iter);
	} // if max_iter > 0

	double fval = the_minimizer->f;
	gsl_multimin_fdfminimizer_free(the_minimizer);
	return fval;
}


/* findMinEuclidianDist - finds the minimum (and max) Euclidian distance between all paris of points
	in the design matrix X (excluding identical design points). After returning,
	m1 points to the minimum distance and m2 points to the maximum distance
*/
void findMinEuclidianDist(gsl_matrix *X, double *m1, double *m2) {
	int i, j, p;
        *m1 = DBL_MAX;
        *m2 = 0;
	double d, diff;
        for (i = 0; i < X->size1 - 1; i++) {
                for (j = i+1; j < X->size1; j++) {
			d = 0.0;
			for (p = 0; p < X->size2; p++) {
				diff = (X->data[i*X->tda+p] - X->data[j*X->tda+p]);
				d += diff*diff;
			}
                        if (d > 0 && d < *m1) {
                                *m1 = d;
                        }
                        if (d > *m2) {
                                *m2 = d;
                        }
                }
        }
}

/****************************************************************************************************
getUnivariateCorRange
	gsl_matrix *X - pointer to design matrix 
	double *m1 - at end of function call, pointer to value of beta for which max correlation btwn
			any 2 observations is 0.65
	double *m2 - at end of function call, pointer to value of beta for which min correlation btwn
			any 2 observations is 0.3

	note: the purpose of this function is to find initial values of correlation parameters so 
		that var-cov matrix is NOT singular (which could happen if correlation btwn any
		2 observations is too high)
****************************************************************************************************/
void getUnivariateCorRange(gsl_matrix *X, double *m1, double *m2) {
   double xmin, xmax;
   findMinEuclidianDist(X, &xmin, &xmax);
   *m1 = -log(.65) / xmin;    
   *m2 = -log(.3) / xmin;     
}


// simplex_steps must be same size as number of parameters to estimate; we do not check this here
// user can also specify: step_size for simplex, step size for gradient

/***************************************************************************************************************
fitGP - the main function that finds mles for gp parameters; we run 'numSimplexTries' simplexes to find 
		a good initial value, then continue with bfgs method
	gsl_matrix *X - pointer to design matrix that determines correlation structure
	gsl_matrix *Y - pointer to 1-column matrix of observations
	int numParams - the number of columns in X 
	int constantMean - 1 if gp should have a constant mean; otherwise mean function will be (1 X) %*% B
		(linear regression in X with intercept added)
	int numSimplexTries - the number of simplexes to run before moving to bfgs method
	int maxSimplexIterations - maximum number of iterations / simplex
	double break_size - stopping criteria for simplex
	double *simplex_steps - vector of initial step sizes for simplex (of length equal to # parameters)
	int BFGS_max_iter - max # iof iterations for bfgs method
	double BFGS_tol - stopping criteria for bfgs method
	double BFGS_steps - initial step size for bfgs method
	int rng_seed - random # seed
	double *nugget - pointer to initial nugget value or nugget matrix
	int nugget_length - 1 for constant nugget, or greater than 1 for nugget matrix, in which case
		it should be equal to the # of observations
	double min_nugget - minimum value f nugget (added to diagonal of var-cov matrix)
	double *estimates - stores values of mles at end of function call
	int verbose - 0 for no status reporting; 1 otherwise
**************************************************************************************************************/
void fitGP(gsl_matrix *X, gsl_matrix *Y, int numParams, int constantMean, 
	int numSimplexTries, int maxSimplexIterations, double break_size, double *simplex_steps,
	int BFGS_max_iter, double BFGS_tol, double BFGS_h, double BFGS_steps, int rng_seed,
	double *nugget, int nugget_length, double min_nugget, double *estimates, int verbose) {

	if (X->size1 != Y->size1) {
		printerr("ERROR: number of observations does not match number of design points\n");
		abort();
	}
	if (nugget_length > 1 && nugget_length != Y->size1) {
		printerr("ERROR: length of nugget matrix must match number of observations, or use constant nugget\n");
		abort();
	}

	time_t curtime;
	struct tm *loctime;
	curtime = time(NULL);
	loctime = localtime(&curtime);

	gsl_matrix *fX;
	if (constantMean > 0) fX = gsl_matrix_alloc(X->size1, 1);
	else fX = gsl_matrix_alloc(X->size1, X->size2 + 1);
	if (constantMean > 0) {
		gsl_matrix_set_all(fX, 1);
	}
	else {
		gsl_matrix *ones = gsl_matrix_alloc(X->size1, 1);
		gsl_matrix_set_all(ones, 1);
		cbind(ones, X, fX);
	}

	gp_params gp;
	gp.X = X;
	gp.fX = fX;	
	gp.Y = Y;
	gp.num_params = X->size2;  // number of regression parameters for mean
	gp.h = BFGS_h;             // for finite difference approximation
	gp.nugget = NULL;          // pointer to nugget array for non-constant nugget 
	if (nugget_length > 1) gp.nugget = nugget;
	gp.verbose = verbose;
	gp.constantNugget = 0;
	if (nugget_length == 1 && nugget[0] <= 0) gp.constantNugget = 1;

	gp.min_nugget = min_nugget;

	/** the vector of parameters to maximize; correlation parameters + nugget **/
	gsl_vector *v;
	if (nugget_length == 1 && nugget[0] <= 0) v = gsl_vector_alloc(X->size2); // corr params only
	else  v = gsl_vector_alloc(X->size2 + 1);    // corr params + nugget term

	gsl_vector *steps = gsl_vector_alloc(v->size);
	int i;
	for (i = 0; i < v->size; i++) {
		steps->data[i] = simplex_steps[i];
	}	

	/** choose initial values **/
	double sig2 =  gsl_stats_variance(Y->data, 1, Y->size1);   

	if (nugget_length == 1 && nugget[0] > 0) {   // we are adding a constant nugget
		double nug;
		if (nugget_length == 1) nug =  *nugget; 
		else nug = 1.0;
		nug = nug / sig2;	// scale the nugget
		v->data[X->size2] = nug;	
	}
	else if (nugget_length > 1) {
		v->data[X->size2] = 1.0 / sig2;
	}	
		

	double m1 = 0, m2;
	/** correlation btwn 2 closest points is btwn m1 and m2, assuming a constant beta **/
	getUnivariateCorRange(X, &m1, &m2);   

	/************ randomly select initial corr params **/
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, rng_seed);
	int simplex_try;
	
	double fval = 0;
	double best_fval = 999999999999999;
	int best_try = 0;
	gsl_vector *best_v = gsl_vector_alloc(v->size);

	for (simplex_try = 0; simplex_try < numSimplexTries; simplex_try++) { 
	
		if (verbose > 0) printout("running simplex # %d...\n", simplex_try+1);
		
	
		for (i = 0; i < X->size2; i++) {
			v->data[i] = m1 + (m2 - m1) * gsl_rng_uniform(rng); 
		}
	
		// put all initial param values on the log scale
		vector_log_subset(v, X->size2);

		fval = optimizeUsingSimplex(v, &gp, steps, maxSimplexIterations, break_size);
		if (verbose > 0) {
			printout("...simplex #%d complete, loglike = %f\n", simplex_try+1, -fval);
		}
		if (simplex_try == 0 || fval < best_fval) {
			best_fval = fval;
			best_try = simplex_try + 1;
			gsl_vector_memcpy(best_v, v);
		}
	}

		if (verbose > 0) {
			printout("\nusing BFGS method from simplex #%d...\n", best_try);
		}

		fval = optimizeUsingGradient(best_v, &gp, BFGS_steps, BFGS_tol, BFGS_max_iter);


	/***************/
	///// REPORT RESULTS ////

	if (verbose) printout("\nmaximum likelihood estimates found, log like =  %f\n", -fval);   // gradient
	
	// put best estimates on actual scale
	vector_exp_check(best_v);

	// get mle's for beta corr parameters	
	gsl_vector *B = gsl_vector_alloc(numParams);
	for (i = 0; i < numParams; i++) {
		B->data[i] = best_v->data[i];
	}	

	gsl_matrix *corr = gsl_matrix_alloc(X->size1, X->size1);
	createCorrMatrix(X, B, corr);

	if (nugget_length == 1 && nugget[0] > 0) {
		addNugget(corr, best_v->data[numParams]);
	}
	else if (nugget_length > 1) {
		addNuggetMatrix(corr, fabs(best_v->data[numParams]), gp.nugget);
	}
	// add minimum nugget
	addNugget(corr, gp.min_nugget);

	gsl_matrix *corrInv = gsl_matrix_alloc(X->size1,X->size1);
	solve(corr, corrInv);

	
	gsl_matrix *bhat = gsl_matrix_alloc(fX->size2, 1);
	calcBhat(fX, corrInv, Y, bhat);
	gsl_matrix *mu_hat = gsl_matrix_alloc(fX->size1, bhat->size2);
	matrix_multiply(fX, bhat, mu_hat);

	sig2 = calcMLESig2(Y, mu_hat, corrInv);
	gsl_matrix_scale(corr, sig2);

	// free random number generator 
	gsl_rng_free(rng);

	// create solution vector
	int index = 0;
	for (i = 0; i < fX->size2; i++) {
		estimates[index] = bhat->data[i];
		index++;
	}
	for (i = 0; i < B->size; i++) {
		estimates[index] = B->data[i];
		index++;
	}
	estimates[index] = sig2;
	index++;
	// if necessary, add the true (unscaled) value of the nugget to the estimates list
	if (nugget_length > 1 || nugget[0] > 0) {
		estimates[index] = fabs(best_v->data[numParams]) * sig2;
	}

	curtime = time(NULL);
	loctime = localtime(&curtime);

		
}

void fitGPTest(double *in_X, double *in_Y, int *numObs, int *numParams) {
//void fitGPTest(double *in_X, double *in_Y) {

	printout("fitGPTest");
	/**
	gsl_matrix *X = gsl_matrix_alloc(*numObs, *numParams);
	createMatrixByCol(in_X, *numObs, *numParams, X);
	gsl_matrix *Y = gsl_matrix_alloc(*numObs, 1);
	createMatrixByCol(in_Y, *numObs, *numParams, Y);
	printMatrix(X, 0);
	***/
}

/*
	fitGPfromR - the function actually called from R, that calls fitgp. See fitgp for parameter
		descriptions
*/
void fitGPFromR (double *in_X, double *in_Y, int *numObs, int *numParams, int *constantMean, 
	int *numSimplexTries, int *maxSimplexIterations, double *break_size, double *simplex_steps, 
	int *maxBFGSIterations, double *BFGS_tol, double *BFGS_h, double *BFGS_steps, int *rng_seed, 
	double *nugget, int *nugget_length, double *min_nugget, double *estimates, int *verbose) {

	gsl_matrix *X = gsl_matrix_alloc(*numObs, *numParams);
	createMatrixByCol(in_X, *numObs, *numParams, X);
	gsl_matrix *Y = gsl_matrix_alloc(*numObs, 1);
	createMatrixByCol(in_Y, *numObs, 1, Y);

	fitGP(X, Y, *numParams, *constantMean, *numSimplexTries, *maxSimplexIterations, *break_size, simplex_steps,
		*maxBFGSIterations, *BFGS_tol, *BFGS_h, *BFGS_steps, *rng_seed, nugget, *nugget_length, *min_nugget, 
		estimates, *verbose);

	gsl_matrix_free(X);
	gsl_matrix_free(Y);
	
}

#endif


