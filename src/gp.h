
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

/*******************************************************************************
gp.h - finds mles of gaussian process parameters that have closed form
	solutions (i.e, B in the mean matrix X'B, and the overall variance sig2;
	also creates a correlation matrix assuming product exponential 
	structure
*******************************************************************************/

#ifndef __GP__
#define __GP__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>

/*********************************************************************************
calcBhat - finds the mle of the mean regression term X'B of the gaussian process;
	gsl_matrix *X - the design matrix for the regression term X'B
	gsl_matrix *Vinv - the inverse variance-covariance matrix specified up to
			   a multiplicative constant
	gsl_matrix *Y - the observations    
	gsl_matrix *bhat - matrix (with size2 = 1) that stores the mle of B
	
	return value: 1 if there is an error inverting X'X; 0 otherwise

	notes: X is for the mean regression function, not the design matrix for 
		comp model runs
**********************************************************************************/
int calcBhat(gsl_matrix *X, gsl_matrix *Vinv, gsl_matrix *Y, gsl_matrix *bhat) {

	gsl_matrix *Xprime = gsl_matrix_alloc(X->size2, X->size1);
	transpose(X, Xprime);

	gsl_matrix *XprimeVinv = gsl_matrix_alloc(Xprime->size1, Vinv->size2);

	matrix_multiply(Xprime, Vinv, XprimeVinv);

	gsl_matrix *XprimeVinvX = gsl_matrix_alloc(Xprime->size1, X->size2);
	matrix_multiply(XprimeVinv, X, XprimeVinvX);	

	gsl_matrix *ans1 = gsl_matrix_alloc(XprimeVinvX->size1, XprimeVinvX->size2);
	if (solve(XprimeVinvX, ans1) == 1) {
		return 1;
	}

	gsl_matrix *ans2 = gsl_matrix_alloc(XprimeVinvX->size1, Y->size2);
	matrix_multiply(XprimeVinv,Y, ans2);

	matrix_multiply(ans1, ans2, bhat);

	gsl_matrix_free(Xprime);
	gsl_matrix_free(XprimeVinv);
	gsl_matrix_free(XprimeVinvX);
	
	gsl_matrix_free(ans1);
	gsl_matrix_free(ans2);
	return 0;
}


/*********************************************************************************
calcMLEsig2 - finds the mle of the overall variance of the GP;
	gsl_matrix *Y - the observations   
	gsl_matrix *mu - mean matrix with same dimensions as Y 
	gsl_matrix *Vinv - the inverse variance-covariance matrix specified up to
			   a multiplicative constant
	return value: the mle of sig2
********************************************************************************/
double calcMLESig2 (gsl_matrix *Y, gsl_matrix *mu, gsl_matrix *Vinv) {

	gsl_matrix *diff   = gsl_matrix_alloc(Y->size1, 1);
	gsl_matrix *diff_t = gsl_matrix_alloc(1,Y->size1);
	gsl_matrix_memcpy(diff, Y);
	gsl_matrix_sub(diff, mu);
	transpose(diff, diff_t);

	gsl_matrix *diff_tinvV = gsl_matrix_alloc(diff_t->size1, Vinv->size2);
	matrix_multiply(diff_t, Vinv, diff_tinvV);
	
	gsl_matrix *ans = gsl_matrix_alloc(1, 1);
	matrix_multiply(diff_tinvV, diff, ans);

	double sig2 = ans->data[0] / (double)Y->size1;
	gsl_matrix_free(diff);
	gsl_matrix_free(diff_t);
	gsl_matrix_free(diff_tinvV);
	gsl_matrix_free(ans);
	return sig2;
}

/*********************************************************************************
addNugget - adds a constant nugget term to the diagnal of the matrix m 
  	    Note: (this assumes m is a square matrix)
*********************************************************************************/
void addNugget(gsl_matrix *m, double nugget) {
	int i;
	for (i = 0; i < m->size1; i++) {
		m->data[i*m->tda + i] += nugget;			
	}
}


/*********************************************************************************
addNuggetMatrix - adds diag(nugget matrix) the matrix m
*********************************************************************************/
void addNuggetMatrix(gsl_matrix *m, double c, double *nugget_matrix) {
	int i;
	for (i = 0; i < m->size1; i++) {
		m->data[i*m->tda + i] += c * nugget_matrix[i];			
	}
}



/*********************************************************************************
logdmvnorm - calculates and returns the log likelihood of observations x given 
	mean matrix mu and var-cov matrix V; -DBL_MAX is returned if the var-cov
	matrix is singular

	note: matrices x and V will be overwritten
*********************************************************************************/
double logdmvnorm(gsl_vector *x, gsl_vector *mu, gsl_matrix *V) {

	int i;
	gsl_set_error_handler_off();

	i = gsl_linalg_cholesky_decomp(V);      
	gsl_set_error_handler(NULL);

	if (i == 1) { // singular matrix!
		return -DBL_MAX;	
	}

	double logdet = logDetFromCholesky(V);

	gsl_vector *I = gsl_vector_alloc(V->size1); // identity vector (b)
	gsl_vector *inv = gsl_vector_alloc(V->size1); // (x)
	
	gsl_matrix *invV = gsl_matrix_alloc(V->size1,V->size1);
	for (i = 0; i < V->size1; i++) {
		gsl_vector_set_zero(I);
		I->data[i] = 1;
		gsl_linalg_cholesky_solve(V, I, inv);
		assignColumn(invV, inv, i);
	}

	gsl_vector *ans = gsl_vector_alloc(x->size);
	
	gsl_vector_sub(x, mu);  // x now = x - mu	
	xprimeA(x, invV, ans);  // ans = (x-mu)' %*% invV

	double d = dotprod(ans, x);  

	// free memory
	gsl_vector_free(I);
	gsl_vector_free(inv);
	gsl_matrix_free(invV);
	gsl_vector_free(ans);
	return	-(x->size / 2.0) * (M_LN2 + M_LNPI) - 0.5* (logdet + d);
}



/*********************************************************************************
createCorrMatrix - creates the correlation matrix for design X where

	corr(z(xi), z(xj)) = exp[sum(-B(x,p - xj,p')**2)], p = 1,2,..,# of inputs
	
	The correlation matrix is returned in corr, which should be the 
		appropriate size when passed in
*********************************************************************************/
void createCorrMatrix(gsl_matrix *X, gsl_vector *B, gsl_matrix *corr) {

	int i, j, p;
	gsl_matrix_set_all(corr, 0.0);	
	double diff;
	double total;

	for (i = 0; i < X->size1; i++) {
	  for(j = i+1; j < X->size1; j++) {
		diff = total = 0.0;
		for (p = 0; p < B->size; p++) {

			diff = (X->data[i*X->tda + p] - X->data[j*X->tda+p]);
			total += -B->data[p]*diff*diff;
		}
		gsl_matrix_set(corr, i, j, exp(total));
		gsl_matrix_set(corr, j, i, exp(total));
	  }
	}
	for (i = 0; i < X->size1; i++) {
		corr->data[i*corr->tda+i] = 1.0;
	}

}
#endif



