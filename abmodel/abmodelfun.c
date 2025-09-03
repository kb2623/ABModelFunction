#include "abmodelfun.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

const double PI = 3.141592653589793238463;

struct point {
	double x, y, z;
};

typedef struct point Point;

struct proteins {
	char * name;
	char * sequence;
	double optimum;
};

static struct proteins instances[] = {
	{"1BXP", "ABBBBBBABBBAB", -5.6104},
	{"1CB3", "BABBBAABBAAAB", -8.4589},
	{"1BXL", "ABAABBAAAAABBABB", -17.3962},
	{"1EDP", "ABABBAABBBAABBABA", -15.0092},
	{"2ZNF", "ABABBAABBABAABBABA", -18.3402},
	{"1EDN", "ABABBAABBBAABBABABAAB", -21.4703},
	{"2H3S", "AABBAABBBBBABBBABAABBBBBB", -21.1519},
	{"1ARE", "BBBAABAABBABABBBAABBBBBBBBBBB", -25.2800},
	{"2KGU", "ABAABBAABABBABAABAABABABABABAAABBB", -52.7165},
	{"1TZ4", "BABBABBAABBAAABBAABBAABABBBABAABBBBBB", -43.0229},
	{"1TZ5", "AAABAABAABBABABBAABBBBAABBBABAABBABBB", -49.3868},
	{"1AGT", "AAAABABABABABAABAABBAAABBABAABBBABABAB", -65.1990},
	{"1CRN", "BBAAABAAABBBBBAABAAABABAAAABBBAAAAAAAABAAABBAB", -92.9853},
	{"2KAP", "BBAABBABABABABBABABBBBABAABABAABBBBBBABBBAABAAABBABBABBAAAAB", -85.5099},
	{"1HVV", "BAABBABBBBBBAABABBBABBABBABABAAAAABBBABAABBABBBABBAABBABBAABBBBBAABBBBBABBB", -95.4475},
	{"1GK4", "ABABAABABBBBABBBABBABBBBAABAABBBBBAABABBBABBABBBAABBABBBBBAABABAAABABAABBBBAABABBBBA", -106.4190},
	{"1PCH", "ABBBAAABBBAAABABAABAAABBABBBBBABAAABBBBABABBAABAAAAAABBABBABABABABBABBAABAABBBAABBAAABA", -156.5250},
	{"2EWH", "AABABAAAAAAABBBAAAAAABAABAABBAABABAAABBBAAAABABAAABABBAAABAAABAAABAABBAABAAAAABAAABABBBABBAAABAABA", -245.5190},
	{"F13", "ABBABBABABBAB", -6.9961},
	{"F21", "BABABBABABBABBABABBAB", -16.5544},
	{"F34", "ABBABBABABBABBABABBABABBABBABABBAB", -31.3455},
	{"F55", "BABABBABABBABBABABBABABBABBABABBABBABABBABABBABBABABBAB", -51.9030},
	{"F89", "ABBABBABABBABBABABBABABBABBABABBABBABABBABABBABBABABBABABBABBABABBABBABABBABABBABBABABBAB", -81.5297}
};
const int no_instances = 23;

void get_pfo_protein(const char * name, char * sequence, double * optimum) {
	int i;
	for (i = 0; i < no_instances; i++) {
		if (strcmp(name, instances[i].name) == 0) {
			strcpy(sequence, instances[i].sequence);
			*optimum = instances[i].optimum;
			return;
		}
	}
	strcpy(sequence, instances[0].sequence);
	*optimum = instances[0].optimum;
}

ABModel* init_abmodel(char * sequence) {
	ABModel * model = (ABModel *) malloc(sizeof(ABModel));
	unsigned int len = strlen(sequence);
	model->seq_len = len;
	model->dim = (len - 2) * 2;
	for (int i = 0; i < len; i++) {
		char c = toupper(sequence[i]);
		if (c =='A' || c == 'B') model->seq[i] = sequence[i];
		else model->seq[i] = 'A';
	}
	return model;
}

double ** get_abmodel_bounds(ABModel * model) {
	double ** B = (double **) malloc(2 *sizeof(double *));
	B[0] = (double *) malloc(model->dim * sizeof(double));
	B[1] = (double *) malloc(model->dim * sizeof(double));
	for (int i = 0; i < model->dim; i++) {
		B[0][i] = 0;
		B[1][i] = 2 * PI;
	}
	return B;
}

void clean_abmodel_bounds(double ** B) {
	free(B[0]);
	free(B[1]);
	free(B);
}

double evaluate_abmodel(ABModel * m, double * x) {
	double e1 = 0, e2 = 0, d, c, p, x, y, z;
	Point * tertiary = (Point *) malloc((m->seq_len) * sizeof(Point));
	double * cos_theta = (double *) malloc((m->seq_len - 2) * sizeof(double));
	double * theta = x;
	double * beta = x[m->seq_len - 2];
	for (int i = 0; i < m->seq_len - 2; i++) cos_theta[i] = cos(theta[i]);
	// Set first two points that are constant
	tertiary[0].x = 0, tertiary[0].y = 0, tertiary[0].z = 0;
	tertiary[1].x = 0, tertiary[1].y = 1, tertiary[1].z = 0;
	// Set third point based on data
	tertiary[2].x = cos_theta[0], tertiary[2].y = 1 + sin(theta[0]), tertiary[2].z = 0;
	for (int i = 0; i < m->seq_len; i++) {
		tertiary[i].x = tertiary[i - 1].x + cos_theta[i - 2]  * cos(beta[i - 3]);
		tertiary[i].y = tertiary[i - 1].y + sin(theta[i - 2]) * cos(beta[i - 3]);
		tertiary[i].z = tertiary[i - 1].z + sin(beta[i - 3]);
	}
	// Main equation when we have 3D points
	for (int i = 0; i < m->seq_len - 2; i++) {
		e1 += 1 - cos_theta[i];
		for (int j = i + 2; j < m->seq_len; j++) {
			double c = 0;
			if (m->seq[i] + m->seq[j] == 130) c = 1;
			else if (m->seq[i] + m->seq[j] == 131) c = -0.5;
			else c = 0.5;
			double x = tertiary[i].x - tertiary[j].x;
			double y = tertiary[i].y - tertiary[j].y;
			double z = tertiary[i].z - tertiary[j].z;
			double d = x * x + y * y + z * z;
			double p = 1 / (d * d * d);
			e2 += p * p - c * p;
		}
	}
	free(tertiary);
	free(cos_theta);
	return e1 / 4 + 4 * e2;
}

/* Code to free the allocated memory */
void clean_abmodel(ABModel * obj) {
	free(obj->seq);
}