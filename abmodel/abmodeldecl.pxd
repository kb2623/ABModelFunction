#!cpython


cdef extern from "abmodelfun.h":
    cdef struct abmodel:
        unsigned int seq_len
        char * seq
        unsigned int dim

    ctypedef abmodel ABModel


cdef extern from "abmodelfun.h":
    void get_pfo_protein(char *, char *, double *)
    ABModel* init_abmodel(char *)
    double ** get_abmodel_bounds(ABModel *)
    void clean_abmodel_bounds(double **)
    double evaluate_abmodel(ABModel *, double *)
    void clean_abmodel(ABModel *)

