import os
import sys

from setuptools import (
    find_packages,
    Extension,
    setup
)

from Cython.Distutils import build_ext


sourcefiles = [
    'abmodel/abmodelfun.c',
    'abmodel/abmodel.pyx',
]

libs = []
if sys.platform != 'win32':
    libs.append('m')

compiler_args = []
if os.environ.get('DEBUG', '0') == '1':
    if sys.platform == 'win32': compiler_args.extend(['/Zi', '/Od'])
    else: compiler_args.extend(['-g', '-O0'])
else:
    compiler_args.extend(['-std=c2x', '-O3', '-march=native'])

link_args = []
if os.environ.get('DEBUG', '0') == '1' and sys.platform == 'win32':
    link_args.append('/DEBUG')

cec2005real = Extension(
    "abmodel.abmodel",
    sources = sourcefiles,
    include_dirs = ["abmodel"],
    language = "c",
    extra_compile_args = compiler_args,
    extra_link_args= link_args,
    libraries = libs
)


setup(
    packages = find_packages(),
    ext_modules = [cec2005real],
    package_data = {'abmodel': [
        'abmodeldecl.pxd',
        'abmodelfun.h'
    ]},
    tests_require = ['pytest']
)
