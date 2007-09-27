
#include "matrix_vector.h"
#include "gp.h"
#include "fit_gp.h"

#include "print.h"

#ifdef __useR__
#include <R_ext/Utils.h>
#endif

int main() {


	gsl_matrix *x = gsl_matrix_alloc(10,1);
	gsl_matrix *y = gsl_matrix_alloc(10,1);
	double d_y[] = {1.8, 2.2, 3.2, 4.9, 5.29, 6.5, 7.8, 8.4, 9.0, 10.17};
	int i = 0;
	for (i = 0; i < 10; i++) {
		x->data[i*x->tda+0] = i;
		y->data[i*y->tda+0] = d_y[i];
	}

	double *stepSize = malloc(sizeof(double)*1);
	*stepSize = 5.0;

	// estimates for Beta, sig2GP, mu
	/***/
	int numEstimates = 3;
	double *estimates = malloc(sizeof(double)*numEstimates);
	for (i = 0; i < numEstimates; i++) {
		estimates[i] = 0.0;
	} 
	/****/	

	double *nugget = malloc(sizeof(double)*1);
	*nugget = 0;

	//fitGP(X,Y,numParams, constantMean, numSimplexTries, maxSimplexIter, break_size, simplex_steps, 
	//	maxBFGSiter, BFGS_tol, BFGS_h, BFGS_steps, seed, nugget, nugget_length, min_nugget, estimates, verbose)	

	fitGP(x,y,1,1,5,100,1e-10, stepSize, 
			500, .01, 1e-10, .001, 0, nugget, 1, 0, estimates, 1);

	printout("fitGP complete: final estimates (meanReg, Beta, sig2, ...):\n");
	for (i = 0; i < numEstimates; i++) {
		printout("%f  ", estimates[i]);
	}
	printout("\n\n");

	return 1;
}



