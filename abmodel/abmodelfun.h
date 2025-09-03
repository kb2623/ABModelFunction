#ifndef ABMODELFUN_H
#define ABMODELFUN_H

struct abmodel {
    unsigned int seq_len;   /* Length of AB sequence */
    char * seq;             /* AB sequence represented with letters A and B */
    unsigned int dim;       /* Problem dimensionality */
};

typedef struct abmodel ABModel;

void get_pfo_protein(char *, char *, double *);
ABModel* init_abmodel(char *);
double ** get_abmodel_bounds(ABModel *);
void clean_abmodel_bounds(double **);
double evaluate_abmodel(ABModel *, double *);
void clean_abmodel(ABModel *);

#endif
