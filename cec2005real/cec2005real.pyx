#!python
#cython: language_level=2, boundscheck=False
import os
import sys
import pkgutil

import cython

from libc.stdlib cimport malloc, free
from libcpp cimport bool

cdef extern from "cec2005.h":
    bool isBound_cec2005()
    void init_cec2005(int nfun, int dim)
    double eval_cec2005(const double *x, int ndim)
    double eval_cec2005_ld(const long double *x, int ndim)
    void getInfo_cec2005(int fun, char *name, double *min, double *max, double *optime)
    void set_directory(char *dir)

cpdef set_function(int fun, int dim):
    init_cec2005(fun, dim)

cdef int m_fun
cdef int m_dim

def _cec2005_eval_func(double[::1] x):
    cdef int dim
    cdef double fitness
    cdef double * sol

    dim = x.shape[0]

    sol = <double *> malloc(dim * cython.sizeof(double))

    if sol is NULL:
        raise MemoryError()

    for i in xrange(dim):
        sol[i] = x[i]

    fitness = eval_cec2005(sol, dim)
    free(sol)
    return fitness

def file_load(data_dir: str, file_name: str):
    if os.path.exists('%s/%s' % (data_dir, file_name)): return
    data = pkgutil.get_data('cec2005real', 'cdatafiles/%s' % file_name)
    with open('%s/%s' % (data_dir, file_name), 'wb') as f: f.write(data)

cpdef get_num_functions(self):
    return 25

cpdef _get_info(int fun, int dim):
    """
    Return the lower bound of the function
    """
    cdef double optimum
    cdef double minvalue, maxvalue
    cdef char* name = <char *> malloc(300)
    optimum = 0
    getInfo_cec2005(fun, name, &minvalue, &maxvalue, &optimum)

    return {'lower': minvalue, 'upper': maxvalue, 'threshold': 1e-8,
           'best': optimum, 'dimension': dim}

cdef class Function:
    cdef int fun
    cdef int dim
    
    def __init__(self, int fun, int dim):
        self.fun = fun
        self.dim = dim
        os.makedirs('cdatafiles', exits_ok=True)
        cdef bytes dir_name = ('%s/cdatafiles' % os.getcwd()).encode()
        set_directory(dir_name)
        # TODO add copying for input files based on chosen function
        if fun is 1:
            file_load('cdatafiles', 'sphere_func_data.txt')
        elif fun is 2:
            file_load('cdatafiles', 'schwefel_102_data.txt')
        elif fun is 3:
            if dim is 2: file_load('cdatafiles', 'elliptic_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'elliptic_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'elliptic_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'elliptic_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'high_cond_elliptic_rot_data.txt')
        elif fun is 4:
            file_load('cdatafiles', 'schwefel_102_data.txt')
        elif fun is 5:
            file_load('cdatafiles', 'schwefel_206_data.txt')
        elif fun is 6:
            file_load('cdatafiles', 'rosenbrock_func_data.txt')
        elif fun is 7:
            if dim is 2: file_load('cdatafiles', 'griewank_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'griewank_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'griewank_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'griewank_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'griewank_func_data.txt')
        elif fun is 8:
            if dim is 2: file_load('cdatafiles', 'ackley_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'ackley_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'ackley_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'ackley_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'ackley_func_data.txt')
        elif fun is 9:
            file_load('cdatafiles', 'rastrigin_func_data.txt')
        elif fun is 10:
            if dim is 2: file_load('cdatafiles', 'rastrigin_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'rastrigin_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'rastrigin_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'rastrigin_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'rastrigin_func_data.txt')
        elif fun is 11:
            if dim is 2: file_load('cdatafiles', 'weierstrass_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'weierstrass_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'weierstrass_M_D50.txt')
            elif dim is 50: file_load('cdatafiles', 'weierstrass_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'weierstrass_data.txt')
        elif fun is 12:
            file_load('cdatafiles', 'schwefel_213_data.txt')
        elif fun is 13:
            file_load('cdatafiles', 'EF8F2_func_data.txt')
        elif fun is 14:
            if dim is 2: file_load('cdatafiles', 'E_ScafferF6_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'E_ScafferF6_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'E_ScafferF6_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'E_ScafferF6_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'E_ScafferF6_func_data.txt')
        elif fun is 15:
            file_load('cdatafiles', 'hybrid_func1_data.txt')
        elif fun in [16, 17]:
            if dim is 2: file_load('cdatafiles', 'hybrid_func1_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'hybrid_func1_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'hybrid_func1_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'hybrid_func1_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func1_data.txt')
        elif fun in [18, 19, 20]:
            if dim is 2: file_load('cdatafiles', 'hybrid_func2_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'hybrid_func2_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'hybrid_func2_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'hybrid_func2_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func2_data.txt')
        elif fun in [21, 22, 23]:
            if dim is 2: file_load('cdatafiles', 'hybrid_func3_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'hybrid_func3_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'hybrid_func3_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'hybrid_func3_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func3_data.txt')
        elif fun in [24, 25]:
            if dim is 2: file_load('cdatafiles', 'hybrid_func4_M_D2.txt')
            elif dim is 10: file_load('cdatafiles', 'hybrid_func4_M_D10.txt')
            elif dim is 30: file_load('cdatafiles', 'hybrid_func4_M_D30.txt')
            elif dim is 50: file_load('cdatafiles', 'hybrid_func4_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func4_data.txt')
        else:
            raise Exception('Function number not defined!!!')
        init_cec2005(fun, dim)

    cpdef info(self):
        return _get_info(self.fun, self.dim)
    
    cpdef get_dim(self):
        return self.dim

    cpdef get_eval_function(self):
        return _cec2005_eval_func
