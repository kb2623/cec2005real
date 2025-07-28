#!cpython

from libcpp cimport bool


cdef extern from "benchfunctions.h":
    cdef struct cec2005data:
        int nfunc
        int nreal
        double C
        double global_bias
        double global_bias
        double *bias
        double *trans_x
        double *basic_f
        double *temp_x1
        double *temp_x2
        double *temp_x3
        double *temp_x4
        double *weight
        double *sigma
        double *lam
        double *norm_x
        double *norm_f
        double **o
        double **g
        double ***l
        double **Af5
        double *Bf5
        double **Af12
        double **Bf12
        double *alphaf12
    
    ctypedef cec2005data CEC2005data


cdef extern from "cec2005.h":
    CEC2005data * init_cec2005(int, int)
    void getInfo_cec2005(int, char*, double*, double*, double*)
    bool isBound_cec2005()
    double eval_cec2005(double*, CEC2005data*)
    void finish_cec2005(CEC2005data*)

