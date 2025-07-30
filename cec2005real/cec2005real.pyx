#!python
#cython: language_level=2, boundscheck=False
import os
import sys
import pkgutil
import cython

from libc.stdlib cimport malloc, free, srand
from libc.time cimport time

from cec2005decl cimport (
    CEC2005data,
    init_cec2005,
    getInfo_cec2005,
    isBound_cec2005,
    eval_cec2005,
    finish_cec2005
)


cdef void file_load(str data_dir, str file_name):
    if not os.path.exists('%s/%s' % (data_dir, file_name)):
        data = pkgutil.get_data('cec2005real', 'cdatafiles/%s' % file_name)
        with open('%s/%s' % (data_dir, file_name), 'wb') as f: f.write(data)


cdef class Function:
    cdef CEC2005data * fdata

    def __init__(self, int fun, int dim):
        os.makedirs('cdatafiles', exist_ok=True)
        if fun == 1:
            file_load('cdatafiles', 'sphere_func_data.txt')
        elif fun == 2:
            file_load('cdatafiles', 'schwefel_102_data.txt')
        elif fun == 3:
            if dim == 2: file_load('cdatafiles', 'elliptic_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'elliptic_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'elliptic_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'elliptic_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'high_cond_elliptic_rot_data.txt')
        elif fun == 4:
            file_load('cdatafiles', 'schwefel_102_data.txt')
        elif fun == 5:
            file_load('cdatafiles', 'schwefel_206_data.txt')
        elif fun == 6:
            file_load('cdatafiles', 'rosenbrock_func_data.txt')
        elif fun == 7:
            if dim == 2: file_load('cdatafiles', 'griewank_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'griewank_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'griewank_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'griewank_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'griewank_func_data.txt')
        elif fun == 8:
            if dim == 2: file_load('cdatafiles', 'ackley_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'ackley_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'ackley_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'ackley_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'ackley_func_data.txt')
        elif fun == 9:
            file_load('cdatafiles', 'rastrigin_func_data.txt')
        elif fun == 10:
            if dim == 2: file_load('cdatafiles', 'rastrigin_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'rastrigin_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'rastrigin_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'rastrigin_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'rastrigin_func_data.txt')
        elif fun == 11:
            if dim == 2: file_load('cdatafiles', 'weierstrass_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'weierstrass_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'weierstrass_M_D50.txt')
            elif dim == 50: file_load('cdatafiles', 'weierstrass_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'weierstrass_data.txt')
        elif fun == 12:
            file_load('cdatafiles', 'schwefel_213_data.txt')
        elif fun == 13:
            file_load('cdatafiles', 'EF8F2_func_data.txt')
        elif fun == 14:
            if dim == 2: file_load('cdatafiles', 'E_ScafferF6_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'E_ScafferF6_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'E_ScafferF6_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'E_ScafferF6_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'E_ScafferF6_func_data.txt')
        elif fun == 15:
            file_load('cdatafiles', 'hybrid_func1_data.txt')
        elif fun == 16 or fun == 17:
            if dim == 2: file_load('cdatafiles', 'hybrid_func1_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'hybrid_func1_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'hybrid_func1_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'hybrid_func1_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func1_data.txt')
        elif fun == 18 or fun == 19 or fun == 20:
            if dim == 2: file_load('cdatafiles', 'hybrid_func2_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'hybrid_func2_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'hybrid_func2_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'hybrid_func2_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func2_data.txt')
        elif fun == 21 or fun == 22 or fun == 23:
            if dim == 2: file_load('cdatafiles', 'hybrid_func3_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'hybrid_func3_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'hybrid_func3_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'hybrid_func3_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func3_data.txt')
        elif fun == 24 or fun == 25:
            if dim == 2: file_load('cdatafiles', 'hybrid_func4_M_D2.txt')
            elif dim == 10: file_load('cdatafiles', 'hybrid_func4_M_D10.txt')
            elif dim == 30: file_load('cdatafiles', 'hybrid_func4_M_D30.txt')
            elif dim == 50: file_load('cdatafiles', 'hybrid_func4_M_D50.txt')
            else: raise Exception('Undefined dimensionality!!!')
            file_load('cdatafiles', 'hybrid_func4_data.txt')
        else:
            raise Exception('Function number not defined!!!')
        self.fdata = init_cec2005(fun, dim)
    
    cpdef info(self):
        cdef double optimum = 0
        cdef double minvalue, maxvalue
        cdef char* name = <char *> malloc(300)
        getInfo_cec2005(self.fdata.nfunc, name, &minvalue, &maxvalue, &optimum)
        return {
            'lower': minvalue,
            'upper': maxvalue,
            'threshold': 1e-8,
            'best': optimum,
            'dimension': self.fdata.nreal
        }
    
    cpdef set_seed(self, unsigned int seed=0):
        if seed == 0: seed = <unsigned int> time(NULL)
        srand(seed)
    
    cpdef eval(self, double[::1] x):
        # Reserve the array to pass to C
        cdef long double * y = <long double *> malloc(self.fdata.nreal * sizeof(long double))
        if y == NULL: raise MemoryError()
        # Copy the original values
        for i in range(self.fdata.nreal): y[i] = x[i]
        # Calculate the fitness value
        cdef long double fx = eval_cec2005(y, self.fdata)
        # Free the memeory
        free(y)
        # Convert and return the value (python cannot work with long double)
        return float(fx)
        
    def __dealloc__(self):
        finish_cec2005(self.fdata)
