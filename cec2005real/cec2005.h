#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "benchfunctions.h"

CEC2005data* init_cec2005(int, int);
long double calc_benchmark_func(long double*, CEC2005data*);
bool isBound_cec2005();
void getInfo_cec2005(int, char*, double*, double*, double*);
void finish_cec2005(CEC2005data*)

#endif