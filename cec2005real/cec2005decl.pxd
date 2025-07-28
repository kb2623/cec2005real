#!cpython


cdef extern from "benchfunctions.h":
    cdef struct cec2005data:
        int nfunc
        int nreal
        long double C
        long double global_bias
        long double global_bias
        long double *bias
        long double *trans_x
        long double *basic_f
        long double *temp_x1
        long double *temp_x2
        long double *temp_x3
        long double *temp_x4
        long double *weight
        long double *sigma
        long double *lam
        long double *norm_x
        long double *norm_f
        long double **o
        long double **g
        long double ***l
        long double **Af5
        long double *Bf5
        long double **Af12
        long double **Bf12
        long double *alphaf12
    
    ctypedef cec2005data CEC2005data


cdef extern from "cec2005.h":
    CEC2005data * init_cec2005(int, int)
    void getInfo_cec2005(int, char*, double*, double*, double*)
    bool isBound_cec2005()
    long double eval_cec2005(long double*, CEC2005data*)
    void finish_cec2005(CEC2005data*)

