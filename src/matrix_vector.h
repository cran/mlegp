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
matrix.h - provides several matrix and vector functions
***************************************************************************/

#ifndef __MATRIX_VECTOR__
#define __MATRIX_VECTOR__


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdio.h>

#include "print.h"

/************************************************************************
createMatrixByRow -  
	d - array of elements in ROW order
	nrow - # of rows
	ncol - # of columns
	m - pointer to a matrix that will store the new matrix, 
		with size1 = nrow and size2 = ncol
************************************************************************/
void createMatrixByRow(double *d, int nrow, int ncol, gsl_matrix *m) {
	int i, j;
	int count = 0;
	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			m->data[i*m->tda+j] = d[count];
			count++;
		}
	}
}

/************************************************************************
createMatrixByCol -  
	d - array of elements in COLUMN order
	nrow - # of rows
	ncol - # of columns
	m - pointer to a matrix that will store the final matrix, 
		with size1 = nrow and size2 = ncol
************************************************************************/
void createMatrixByCol(double *d, int nrow, int ncol, gsl_matrix *m) {
	int i, j;
	int count = 0;
	for (j = 0; j < ncol; j++) {
		for (i = 0; i < nrow; i++) {
			m->data[i*m->tda+j] = d[count];
			count++;
		}
	}
}


/************************************************************************
assignColumn - overwrites column 'col+1' of matrix 'm' with vector 'v', 
************************************************************************/
void assignColumn(gsl_matrix *m, const gsl_vector *v, int col) {
	int i;
	for (i = 0; i < m->size1; i++) {
		m->data[i*m->tda + col] = v->data[i];
	}
}

/************************************************************************
solveFromCholesky - stores the inverse of a matrix A based on its 
	cholesky decomposition in the matrix ans
	
	gsl_matrix *cholesky - pointer to cholesky decomposition of A
	gsl_matrix *ans - pointer to matrix that stores the inverse
************************************************************************/
void solveFromCholesky(gsl_matrix *cholesky, gsl_matrix *ans) {
	// we solve Ax = b for all rows of A 
	gsl_vector *I = gsl_vector_alloc(cholesky->size1); // identity vector (b)
	gsl_vector *inv = gsl_vector_alloc(cholesky->size1); // (x)
	int i;
	for (i = 0; i < cholesky->size1; i++) {
		gsl_vector_set_zero(I);
		I->data[i] = 1;
		gsl_linalg_cholesky_solve(cholesky, I, inv);
		assignColumn(ans, inv, i);
	}
	gsl_vector_free(I);
	gsl_vector_free(inv);
}

/*************************************************************************
solve - stores the inverse of a matrix V in the matrix ans
	gsl_matrix *V - pointer to matrix to invert
	gsl_matrix *ans - will store inverse of V; must be same size as V

	returns 1 if V cannot be inverted; (probably) 0 otherwise
*************************************************************************/
int solve(gsl_matrix *V, gsl_matrix *ans) {

	gsl_set_error_handler_off();
	gsl_matrix *V_copy = gsl_matrix_alloc(V->size1, V->size2);
	gsl_matrix_memcpy(V_copy, V);
	int i = gsl_linalg_cholesky_decomp(V_copy);
	gsl_set_error_handler(NULL);
	if (i == 1) return i;
	solveFromCholesky(V_copy, ans);
	gsl_matrix_free(V_copy);
	return i;

}

// converts matrix to vector, row by row
void toVector(gsl_matrix *m, gsl_vector *v) {
	int i, j, count;
	count = 0;
	for (i = 0; i < m->size1; i++) {
		for (j = 0; j < m->size2; j++) {
			v->data[count] = m->data[i*m->tda + j];
			count++;
		}
	}
}

/****************************************************************************************
detFromCholesky - returns the determinant of matrix A given its cholesky decomposition
****************************************************************************************/
double detFromCholesky(gsl_matrix *ch) {
	double ans = 1.0;
	int i = 0;
	for (i = 0; i < ch->size1; i++) {
		ans *= ch->data[i*ch->tda + i];
	}
	return ans*ans;
}

/*********************************************************************************************
logDetFromCholesky - returns the log determinant of matrix A given its cholesky decomposition
*********************************************************************************************/
double logDetFromCholesky(gsl_matrix *ch) {
	double ans = 0.0;
	int i = 0;
	for (i = 0; i < ch->size1; i++) {
		ans += log(ch->data[i*ch->tda + i]);
	}
	return 2.0*ans;
}


/*********************************************************************************************
matrix_multiply - stores A%*%B in the matrix ans
*********************************************************************************************/
void matrix_multiply(gsl_matrix *A, gsl_matrix *B, gsl_matrix *ans) {

	if (A->size2 != B->size1) {
		printerr("ERROR: matrices are not compatible for matrix multiplication!\n");
		abort();
	}
	int i, j, p;
	double total;
	for (i = 0; i < A->size1; i++) {
		for (j = 0; j < B->size2; j++) {
			total = 0.0000;
			for (p = 0; p < A->size2; p++) {
				total += A->data[i*A->tda+p] * B->data[p*B->tda + j];
			}
			ans->data[i*ans->tda + j] = total;
		}
	}
}

/***********************************************************************
cbind - implementation of cbind from R; stores the matrix {X1 X2} in A
***********************************************************************/
void cbind(gsl_matrix *X1, gsl_matrix *X2, gsl_matrix *A) {
	int i;
	int j;
	int end = X1->size2 + X2->size2;

	for (i = 0; i < A->size1; i++) {
		for (j = 0; j < X1->size2; j++) {
			A->data[i*A->tda+j] = X1->data[i*X1->tda+j];
		}
		for (j = X1->size2; j< end; j++) {
			A->data[i*A->tda+j] = X2->data[i*X2->tda+j-X1->size2];
		}
	}

}


/********************************************************
transpose - stores transpose of matrix A in Aprime 
********************************************************/
void transpose(gsl_matrix *A, gsl_matrix *Aprime) {
	int i, j;
	for (i = 0; i < A->size1; i++) {
		for (j = 0; j < A->size2; j++) {
			Aprime->data[j*Aprime->tda + i] = A->data[i*A->tda+j];			
		}
	}
}

/*****************************************************************
printMatrix - prints line # 'line' of matrix m, or entire matrix
		if 'line' = 0
*****************************************************************/
void printMatrix(gsl_matrix *m, int line) {
  int i = 0;
  int j = 0;
  if (line != 0) { 
	for (j = 0; j < m->size2; j++) {
		printout("%f ", m->data[(line-1)*m->tda + j]);
	}
	printout("\n");
	return;
  }
   for (i = 0; i < m->size1; i++) {
	for (j = 0; j < m->size2; j++) {
		printout("%f ", m->data[i*m->tda + j]);
	}
	printout("\n");
  }

}


/*****************************************************************
printVector - prints all elements of vector v to 1 row of output
*****************************************************************/
void printVector(gsl_vector *v) {
	int i;
	for (i = 0; i < v->size; i++) {
		printout("%f ", v->data[i]);
	}
	printout("\n");
}


/***********************************************************************
xprimeA - stores the matrix x'A in the vector ans
	  gsl_vector *x - pointer to vector of length n
	  gsl_matrix *A - pointer to matrix of length n,m
	  gsl_vector *ans - stores the answer; must be length m
***********************************************************************/
void xprimeA(gsl_vector *x, gsl_matrix *A, gsl_vector *ans) {
	int i,j;
	double total = 0.0;
	for (j = 0; j < A->size2; j++) {
		total = 0.0;
		for (i = 0; i < A->size1; i++) {
			 total += A->data[i*A->tda+j]*x->data[i];
		}
		ans->data[j] = total;
	}
}


/**************************************************************************
dotprod - calculates and returns the dot product of vectors v1 and v2; 
	v1 will be overwritten
**************************************************************************/
double dotprod(gsl_vector *v1, gsl_vector *v2) {
	int i;
	double total = 0.0;
	gsl_vector_mul(v1, v2);
	for (i = 0; i < v1->size; i++) {
		total += v1->data[i];
	}
	return total;
}

#endif



