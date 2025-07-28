#ifndef UTIL_H
#define UTIL_H

/* Global Constants */
#define INF MAXDOUBLE

double maximum (double, double);
double minimum (double, double);
double modulus (double*, int);
double dot (double*, double*, int);
double mean (double*, int);

void transform(double*, int, int, double*, double*, double*, double*, double*, double**, double**, double***);
void transform_norm(int, int, double*, double*, double*, double*, double**, double***);
void calc_weight(double*, int, int, double*, double*, double**);

#endif