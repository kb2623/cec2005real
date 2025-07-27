#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/* Global Constants */
#define EPS 1.0e-10
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

long double calc_ackley (long double*, int);
long double calc_rastrigin (long double*, int);
long double calc_weierstrass (long double*, int);
long double calc_griewank (long double*, int);
long double calc_sphere (long double*, int);
long double calc_schwefel (long double*, int);
long double calc_rosenbrock (long double*, int);
long double nc_rastrigin (long double*, int, long double*);
long double nc_schaffer (long double, long double);
long double ExpandedF6 (long double, long double)

#endif