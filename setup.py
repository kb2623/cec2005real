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
        '*.h'
        'cec2005real/def1.c',
        'cec2005real/def2.c',
        'cec2005real/def3.c',
        'cec2005real/def4.c',
        'cec2005real/opr_funcs.c',
        'cec2005real/rand.c'
    ]},
    tests_require = ['pytest']
)
