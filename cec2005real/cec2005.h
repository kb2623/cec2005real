#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "benchfunctions.h"

CEC2005data* init_cec2005(int, int);
void getInfo_cec2005(int, char*, double*, double*, double*);
bool isBound_cec2005();
long double eval_cec2005(long double*, CEC2005data*);
void finish_cec2005(CEC2005data*)

#endif