#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/* Global Constants */
#define EPS 1.0e-10
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

double calc_ackley (double*, int);
double calc_rastrigin (double*, int);
double calc_weierstrass (double*, int);
double calc_griewank (double*, int);
double calc_sphere (double*, int);
double calc_schwefel (double*, int);
double calc_rosenbrock (double*, int);
double nc_rastrigin (double*, int, double*);
double nc_schaffer (double, double);
double ExpandedF6 (double, double);

#endif