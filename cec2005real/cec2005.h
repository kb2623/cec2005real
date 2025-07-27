#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "benchfunctions.h"

bool isBound_cec2005();
void init_cec2005(int, int);
double eval_cec2005(const double*, int);
double eval_cec2005_ld(const long double*, int);
void getInfo_cec2005(int, char*, double*, double*, double*);
void set_directory(char*);

#endif