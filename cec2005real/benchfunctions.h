#ifndef BENCHFUNCTIONS_H
#define BENCHFUNCTIONS_H

struct cec2005data {
    int nfunc;              // number of the benchmark function to optimize
    int nreal;              // number of components/variables for functions
    double C;
    double global_bias;
    double *bias;      // functions error fix
    double *trans_x;   // Transformation vector of vector x
    double *basic_f;   // Vector of objective values
    double *temp_x1;   // Temporary vector
    double *temp_x2;   // Temporary vector
    double *temp_x3;   // Temporary vector
    double *temp_x4;   // Temporary vector
    double *weight;    // Vector of weights for functions
    double *sigma;
    double *lam;
    double *norm_x;
    double *norm_f;
    double **o;
    double **g;
    double ***l;
    double **Af5;      // First temporary vector for f5
    double *Bf5;       // Second temporary vector for f5
    double **Af12;     // First temporary vector for f12
    double **Bf12;     // Second temporary vector for f12
    double *alphaf12;  // Third temporary vector for f12
};

typedef struct cec2005data CEC2005data;

double calc_benchmark_func_f1(double*, CEC2005data*);
double calc_benchmark_func_f2(double*, CEC2005data*);
double calc_benchmark_func_f3(double*, CEC2005data*);
double calc_benchmark_func_f4(double*, CEC2005data*);
double calc_benchmark_func_f5(double*, CEC2005data*);
double calc_benchmark_func_f6(double*, CEC2005data*);
double calc_benchmark_func_f7(double*, CEC2005data*);
double calc_benchmark_func_f8(double*, CEC2005data*);
double calc_benchmark_func_f9(double*, CEC2005data*);
double calc_benchmark_func_f10(double*, CEC2005data*);
double calc_benchmark_func_f11(double*, CEC2005data*);
double calc_benchmark_func_f12(double*, CEC2005data*);
double calc_benchmark_func_f13(double*, CEC2005data*);
double calc_benchmark_func_f14(double*, CEC2005data*);
double calc_benchmark_func_f15(double*, CEC2005data*);
double calc_benchmark_func_f16(double*, CEC2005data*);
double calc_benchmark_func_f17(double*, CEC2005data*);
double calc_benchmark_func_f18(double*, CEC2005data*);
double calc_benchmark_func_f19(double*, CEC2005data*);
double calc_benchmark_func_f20(double*, CEC2005data*);
double calc_benchmark_func_f21(double*, CEC2005data*);
double calc_benchmark_func_f22(double*, CEC2005data*);
double calc_benchmark_func_f23(double*, CEC2005data*);
double calc_benchmark_func_f24(double*, CEC2005data*);
double calc_benchmark_func_f25(double*, CEC2005data*);


void calc_benchmark_norm_f15(CEC2005data*);
void calc_benchmark_norm_f16(CEC2005data*);
void calc_benchmark_norm_f17(CEC2005data*);
void calc_benchmark_norm_f18(CEC2005data*);
void calc_benchmark_norm_f19(CEC2005data*);
void calc_benchmark_norm_f20(CEC2005data*);
void calc_benchmark_norm_f21(CEC2005data*);
void calc_benchmark_norm_f22(CEC2005data*);
void calc_benchmark_norm_f23(CEC2005data*);
void calc_benchmark_norm_f24(CEC2005data*);
void calc_benchmark_norm_f25(CEC2005data*);

#endif