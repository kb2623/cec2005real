#include "util.h"

#include <math.h>

/* Function to return the maximum of two variables */
long double maximum (long double a, long double b) {
	if (a>b) {
		return(a);
	}
	return (b);
}

/* Function to return the minimum of two variables */
long double minimum (long double a, long double b) {
	if (a<b) {
		return (a);
	}
	return (b);
}

/* Function to return the modulus of a vector */
long double modulus (long double *x, int n) {
	int i;
	long double res;
	res = 0.0;
	for (i=0; i<n; i++) {
		res += x[i]*x[i];
	}
	return sqrt(res);
}

/* Function to return the dot product of two vecors */
long double dot (long double *a, long double *b, int n) {
	int i;
	long double res;
	res = 0.0;
	for (i=0; i<n; i++) {
		res += a[i]*b[i];
	}
	return (res);
}

/* Function to return the mean of n variables */
long double mean (long double *x, int n) {
	int i;
	long double res;
	res = 0.0;
	for (i=0; i<n; i++) {
		res += x[i];
	}
	return (res / (long double)n);
}

/* Code to transform a variable vector based on function index 'count' */
void transform(long double *x, int count, int nreal, long double *temp_x1, long double *temp_x2, long double *temp_x3, long double *trans_x, long double *lambda, long double **o, long double **g, long double ***l) {
	int i, j;
	for (i = 0; i < nreal; i++) {
		temp_x1[i] = x[i] - o[count][i];
	}
	for (i = 0; i < nreal; i++) {
		temp_x2[i] = temp_x1[i] / lambda[count];
	}
	for (j = 0; j < nreal; j++) {
		temp_x3[j] = 0.0;
		for (i = 0; i < nreal; i++) {
			temp_x3[j] += g[i][j] * temp_x2[i];
		}
	}
	for (j = 0; j < nreal; j++) {
		trans_x[j] = 0.0;
		for (i = 0; i < nreal; i++) {
			trans_x[j] += l[count][i][j] * temp_x3[i];
		}
	}
	return;
}

/* Code to transform a vector (with elements 5.0) based on function index 'count' */
void transform_norm(int count, int nreal, long double *temp_x2, long double *temp_x3, long double *trans_x, long double *lambda, long double **g, long double ***l) {
	int i, j;
	for (i = 0; i < nreal; i++) {
		temp_x2[i] = 5.0 / lambda[count];
	}
	for (j = 0; j < nreal; j++) {
		temp_x3[j] = 0.0;
		for (i = 0; i < nreal; i++) {
			temp_x3[j] += g[i][j] * temp_x2[i];
		}
	}
	for (j = 0; j < nreal; j++) {
		trans_x[j] = 0.0;
		for (i = 0; i < nreal; i++) {
			trans_x[j] += l[count][i][j] * temp_x3[i];
		}
	}
	return;
}

/* Code to compute the weights for a variable vector */
void calc_weight(long double *x, int nfunc, int nreal, long double *weight, long double *sigma, long double **o) {
	int i, j;
	long double sum;
	long double max;
	max = -INF;
	for (i = 0; i < nfunc; i++) {
		sum = 0.0;
		for (j = 0; j < nreal; j++) {
			sum += (x[j] - o[i][j]) * (x[j] - o[i][j]);
		}
		weight[i] = exp(-(sum) / (2.0 * nreal * sigma[i] * sigma[i]));
		max = maximum(max, weight[i]);
	}
	sum = 0.0;
	for (i = 0; i < nfunc; i++) {
		if (weight[i] != max) {
			weight[i] *= (1.0 - pow(max, 10.0));
		}
		sum += weight[i];
	}
	if (sum == 0.0) {
		for (i = 0; i < nfunc; i++) {
			weight[i] = 1.0 / (long double)nfunc;
		}
	} else {
		for (i = 0; i < nfunc; i++) {
			weight[i] /= sum;
		}
	}
	return;
}