#!python
#cython: language_level=2, boundscheck=False
import os
import sys
import pkgutil
import cython

from libc.stdlib cimport malloc, free

from abmodeldecl cimport (
    ABModel,
    get_pfo_protein,
    init_abmodel,
    get_abmodel_bounds,
    clean_abmodel_bounds,
    evaluate_abmodel,
    clean_abmodel
)


def pfo_instance(name: str) -> tuple[str, float]:
    cdef bytes cname = name.encode('utf-8')
    cdef char * seq = <char *> malloc(300)
    cdef double opt
    get_pfo_protein(cname, seq, &opt)
    return seq.decode('utf-8'), opt


cdef class Model:
    cdef ABModel * model

    def __init__(self, str seq):
        cdef bytes b = seq.encode('utf-8')
        cdef char * cseq = b
        self.model = init_abmodel(cseq)

    def get_bounds(self) -> tuple[list, list]:
        cdef double ** B = get_abmodel_bounds(self.model)
        Bl, Bu = [], []
        for i in range(self.model.dim): Bl.append(B[0][i])
        for i in range(self.model.dim): Bu.append(B[1][i])
        clean_abmodel_bounds(B)
        return Bl, Bu
    
    cpdef float eval(self, double[::1] x):
        # Reserve the array to pass to C
        cdef double * y = <double *> malloc(self.model.dim * sizeof(double))
        if y == NULL: raise MemoryError()
        # Copy the original values
        for i in range(self.model.dim): y[i] = x[i]
        # Calculate the fitness value
        cdef double fx = evaluate_abmodel(self.model, y)
        # Free the memeory
        free(y)
        # Return the value
        return fx
        
    def __dealloc__(self) -> None:
        clean_abmodel(self.model)
