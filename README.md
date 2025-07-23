# Introduction

This is a Python wrapping using the C++ Implementation of the test suite for the Special Session on Large Scale Global Optimization at 2005 IEEE Congress on Evolutionary Computation.

## Note

If you are to use any part of this code, please cite the following publications:

P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y.-P. Chen, A. Auger and S. Tiwari, "Problem Definitions and Evaluation Criteria for the CEC 2005 Special Session on Real-Parameter Optimization", Technical Report, Nanyang Technological University, Singapore, May 2005 AND KanGAL Report #2005005, IIT Kanpur, India.

http://www3.ntu.edu.sg/home/EPNSugan/index_files/CEC-05/CEC05.htm   

## Requirements

* GNU gcc
* Python
* Cython

## Testing Environment

* buntu noble/sid
* gcc 13.
* Python 3.13

## Quickstart

The package is very simple to use. There is a class Function with two functions:

* Give information for each function: their optimum, their dimensionality, the domain search, and the expected threshold to achieve the optima.
* Give a fitness function to evaluate solutions. It expect that these solutions are numpy arrays (vectors) but it can also work with normal arrays.

These two functionalities are done with two methods in Benchmark class:

- **get_num_functions()**

  Return the number of functions in the benchmarks (15)

- **get_info()**

  Return an array with the following information, where /function_id/ is the identifier of the function, a int value between 1 and 15.

    - lower, upper
        *lower* and *upper* boundaries of the domain search. 

    - best
        Optimum to achieve, it is always zero, thus it can be ignored.

    - threshold
        Threshold to obtain, it is always zero, thus it can also be ignored.

    - dimension
        Dimension for the function, it is always 1000.

    It can be noticed that several data are the same for all functions. It is made for maintaining the 
    same interface to other cec20xx competitions.

- **get_eval_function()**

  It returns the fitness function to evaluate the solutions.

## Examples of use

### Obtain information about one function

```python
>>> from cec2005real.cec2005 import Function
>>> fbench = Function(1, 10)
>>> fbench.get_info()
{'best': 0.0,
 'dimension': 1000,
 'lower': -100.0,
 'threshold': 0,
 'upper': 100.0}
```

### Create random solution for the search

```python
>>> from numpy.random import rand
>>> info = fbench.get_info()
>>> dim = info['dimension']
>>> sol = info['lower']+rand(dim)*(info['upper']-info['lower'])
```

### Evaluate a solution

```python
>>> fun_fitness = fbench.get_eval_function()
>>> fun_fitness(sol)
464006824710.75995
```
