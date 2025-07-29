#ifndef UTIL_H
#define UTIL_H

#include <float.h>

/* Global Constants */
#define INF DBL_MAX

long double maximum (long double, long double);
long double minimum (long double, long double);
long double modulus (long double*, int);
long double dot (long double*, long double*, int);
long double mean (long double*, int);

void transform(long double*, int, int, long double*, long double*, long double*, long double*, long double*, long double**, long double**, long double***);
void transform_norm(int, int, long double*, long double*, long double*, long double*, long double**, long double***);
void calc_weight(long double*, int, int, long double*, long double*, long double**);

#endif
