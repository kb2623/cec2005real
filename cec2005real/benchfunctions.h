#ifndef BENCHFUNCTIONS_H
#define BENCHFUNCTIONS_H

struct cec2005fdata {
    int nfunc; // number of functions in a benchmark
    int nreal; // number of components/variables for functions
    long double C;
    long double global_bias;
    long double *bias; // functions error fix
    long double *trans_x;
    long double *basic_f;
    long double *temp_x1;
    long double *temp_x2;
    long double *temp_x3;
    long double *temp_x4;
    long double *weight;
    long double *sigma;
    long double *lambda;
    long double *norm_x;
    long double *norm_f;
    long double **o;
    long double **g;
    long double ***l;
    long double **Af5;
    long double *Bf5;
    long double **Af12;
    long double **Bf12;
    long double *alphaf12;
};

typedef struct cec2005fdata CEC2005data;

long double calc_benchmark_func_f1(long double*, CEC2005data*);
long double calc_benchmark_func_f2(long double*, CEC2005data*);
long double calc_benchmark_func_f3(long double*, CEC2005data*);
long double calc_benchmark_func_f4(long double*, CEC2005data*);
long double calc_benchmark_func_f5(long double*, CEC2005data*);
long double calc_benchmark_func_f6(long double*, CEC2005data*);
long double calc_benchmark_func_f7(long double*, CEC2005data*);
long double calc_benchmark_func_f8(long double*, CEC2005data*);
long double calc_benchmark_func_f9(long double*, CEC2005data*);
long double calc_benchmark_func_f10(long double*, CEC2005data*);
long double calc_benchmark_func_f11(long double*, CEC2005data*);
long double calc_benchmark_func_f12(long double*, CEC2005data*);
long double calc_benchmark_func_f13(long double*, CEC2005data*);
long double calc_benchmark_func_f14(long double*, CEC2005data*);
long double calc_benchmark_func_f15(long double*, CEC2005data*);
long double calc_benchmark_func_f16(long double*, CEC2005data*);
long double calc_benchmark_func_f17(long double*, CEC2005data*);
long double calc_benchmark_func_f18(long double*, CEC2005data*);
long double calc_benchmark_func_f19(long double*, CEC2005data*);
long double calc_benchmark_func_f20(long double*, CEC2005data*);
long double calc_benchmark_func_f21(long double*, CEC2005data*);
long double calc_benchmark_func_f22(long double*, CEC2005data*);
long double calc_benchmark_func_f23(long double*, CEC2005data*);
long double calc_benchmark_func_f24(long double*, CEC2005data*);
long double calc_benchmark_func_f25(long double*, CEC2005data*);

#endif