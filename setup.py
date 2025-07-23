import sys

from setuptools import (
    find_packages,
    Extension,
    setup
)

from Cython.Distutils import build_ext

sourcefiles = [
    'cec2005real/cec2005real.pyx',
    'cec2005real/cec2005.c'
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
    name = 'cec2005real',
    packages = find_packages(),
    ext_modules = [cec2005real],
    package_data = {'cec2005real': [
        'cdatafiles/*.txt', 
        'cec2005.h',
        'global.h',
        'interface.h',
        'rand.h',
        'sprng.h',
        'sub.h',
        'def1.c',
        'def2.c',
        'def3.c',
        'def4.c',
        'opr_funcs.c',
        'rand.c'
    ]},
    tests_require = ['pytest']
)
