#include "mrandom.h"

#include <stdlib.h>
#include <math.h>

/**
 * Implementation of normal random generator.
 *
 * To set the seed use srandom() with parameter that represents the seed
 */
double randomnormaldeviate() {
	double x = (double)rand() / RAND_MAX, y = (double)rand() / RAND_MAX;
	return sqrt(-2 * log(x)) * cos(2 * M_PI * y);
}

/*
 * Normal random numbers generator - Marsaglia algorithm.
 */
double *generate(int n) {
	int i;
	int m = n + n % 2;
	double *values = (double *)calloc(m, sizeof(double));
	double average, deviation;
	if (values) {
		for (i = 0; i < m; i += 2) {
			double x, y, rsq, f;
			do {
				x = 2.0 * rand() / (double)RAND_MAX - 1.0;
				y = 2.0 * rand() / (double)RAND_MAX - 1.0;
				rsq = x * x + y * y;
			} while (rsq >= 1. || rsq == 0.);
			f = sqrt(-2.0 * log(rsq) / rsq);
			values[i] = x * f;
			values[i + 1] = y * f;
		}
	}
	return values;
}

/*
void msrand(unsigned int seed) {
	holdrand = (long)seed;
}

int mrand(void) {
	return (((holdrand = holdrand * 214013L + 2531011L) >> 16) & 0x7fff);
}
*/