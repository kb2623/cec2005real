#!python
#cython: language_level=2, boundscheck=False
import os
import sys
import pkgutil
import cython

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
from libcpp cimport bool

from cec2005decl cimport (
    CEC2005data,
    init_cec2005,
    getInfo_cec2005,
    isBound_cec2005,
    eval_cec2005,
    finish_cec2005
)

ctypedef np.longdouble_t np_ldouble


def file_load(data_dir: str, file_name: str):
    if os.path.exists('%s/%s' % (data_dir, file_name)): return
    data = pkgutil.get_data('cec2005real', 'cdatafiles/%s' % file_name)
    with open('%s/%s' % (data_dir, file_name), 'wb') as f: f.write(data)


cpdef _get_info(int fun, int dim):
    """
    Return the lower bound of the function
    """
    cdef double optimum
    cdef double minvalue, maxvalue
    cdef char* name = <char *> malloc(300)
    optimum = 0
    getInfo_cec2005(fun, name, &minvalue, &maxvalue, &optimum)
    return {
        'lower': minvalue,
        'upper': maxvalue,
        'threshold': 1e-8,
        'best': optimum,
        'dimension': dim
    }


cdef class Function:
    cdef CEC2005data * fdata
    cdef bint initialized

    def __cinit__(self, int fun, int dim):
        self.fdata = init_cec2005(fun, dim)
        self.initialized = self.fdata is not NULL
    
    def __init__(self, int fun, int dim):
        os.makedirs('cdatafiles', exist_ok=True)
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

    cpdef info(self):
        if not self.initialized: raise ValueError('Data not initialized')
        return _get_info(self.fdata.nfun, self.fdata.nreal)
    
    cpdef eval(self, x):
        if not self.initialized: raise ValueError('Data not initialized')
        x_arr = np.asarray(x, dtype=np.longdouble)
        # Check dimensionality
        if x_arr.ndim != 1: raise ValueError("Input must be a 1D array or list.")
        if x_arr.shape[0] != self.cdata.nreal: raise ValueError(f"Input vector must have dimension {self.cdata.nreal}.")
        x_arr = np.ascontiguousarray(x_arr, dtype=np.longdouble)
        cdef np.ndarray[np_ldouble, ndim=1, mode="c"] x_np = x_arr
        cdef long double* x_ptr = <long double*> x_np.data
        cdef long double result = eval_cec2005(x_ptr, self.cdata)
        return float(result)

    def __dealloc__(self):
        if self.initialized: finish_cec2005(self.fdata)
