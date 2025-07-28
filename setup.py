import sys

from setuptools import (
    find_packages,
    Extension,
    setup
)

from Cython.Distutils import build_ext

sourcefiles = [
    'cec2005real/util.c',
    'cec2005real/mrandom.c',
    'cec2005real/functions.c',
    'cec2005real/benchfunctions.c',
    'cec2005real/cec2005.c',
    'cec2005real/cec2005real.pyx',
]


cec2005real = Extension(
    "cec2005real.cec2005real",
    sources = sourcefiles,
    include_dirs = ["cec2005real"],
    language = "c",
    extra_compile_args = [
        "-std=c2x",
        "-O3"
    ],
    libraries = ["m"]
)


setup(
    packages = find_packages(),
    ext_modules = [cec2005real],
    package_data = {'cec2005real': [
        'cdatafiles/*.txt', 
        'cec2005decl.pxd',
        'util.h',
        'mrandom.h',
        'functions.h',
        'benchfunctions.h',
        'cec2005.h'
    ]},
    tests_require = ['pytest']
)
