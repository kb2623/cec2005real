#define _ALL

#include "cec2005.h"

#include <assert.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* Code to allocate memory to global variables being used in evaluation of functions */
CEC2005data* allocate_memory(int nreal, int nfunc) {
	int i, j, k;
	CEC2005data *obj = (CEC2005data*)malloc(sizeof(CEC2005data))
	obj->nreal = nreal;
	obj->nfunc = nfunc;
	obj->norm_x = (long double *)malloc(nreal * sizeof(long double));
	obj->norm_f = (long double *)malloc(nfunc * sizeof(long double));
	obj->trans_x = (long double *)malloc(nreal * sizeof(long double));
	obj->basic_f = (long double *)malloc(nfunc * sizeof(long double));
	obj->temp_x1 = (long double *)malloc(nreal * sizeof(long double));
	obj->temp_x2 = (long double *)malloc(nreal * sizeof(long double));
	obj->temp_x3 = (long double *)malloc(nreal * sizeof(long double));
	obj->temp_x4 = (long double *)malloc(nreal * sizeof(long double));
	obj->weight = (long double *)malloc(nfunc * sizeof(long double));
	obj->sigma = (long double *)malloc(nfunc * sizeof(long double));
	obj->lam = (long double *)malloc(nfunc * sizeof(long double));
	obj->bias = (long double *)malloc(nfunc * sizeof(long double));
	obj->o = (long double **)malloc(nfunc * sizeof(long double));
	obj->l = (long double ***)malloc(nfunc * sizeof(long double));
	obj->g = (long double **)malloc(nreal * sizeof(long double));
	for (i = 0; i < nfunc; i++) {
		obj->o[i] = (long double *)malloc(nreal * sizeof(long double));
		obj->l[i] = (long double **)malloc(nreal * sizeof(long double));
	}
	for (i = 0; i < nreal; i++) {
		obj->g[i] = (long double *)malloc(nreal * sizeof(long double));
	}
	for (i = 0; i < nfunc; i++) for (j = 0; j < nreal; j++) {
		obj->l[i][j] = (long double *)malloc(nreal * sizeof(long double));

	}
	/* Do some trivial (common) initialization here itself */
	obj->C = 2000.0;
	for (i = 0; i < nreal; i++) {
		obj->norm_x[i] = 5.0;
		obj->trans_x[i] = 0.0;
		obj->temp_x1[i] = 0.0;
		obj->temp_x2[i] = 0.0;
		obj->temp_x3[i] = 0.0;
		obj->temp_x4[i] = 0.0;
		for (j = 0; j < nreal; j++) {
			if (i == j) {
				obj->g[i][j] = 1.0;
			} else {
				obj->g[i][j] = 0.0;
			}
		}
	}
	for (i = 0; i < nfunc; i++) {
		obj->basic_f[i] = 0.0;
		obj->norm_f[i] = 0.0;
		obj->weight[i] = 1.0 / (long double)nfunc;
		obj->sigma[i] = 1.0;
		obj->lam[i] = 1.0;
		obj->bias[i] = 100.0 * (long double)i;
		for (j = 0; j < nreal; j++) {
			o[i][j] = 0.0;
			for (k = 0; k < nreal; k++) {
				if (j == k) {
					obj->l[i][j][k] = 1.0;
				} else {
					obj->l[i][j][k] = 0.0;
				}
			}
		}
	}
	return obj;
}

/* Code to free the allocated memory */
void free_memory(CEC2005data *obj) {
	int i, j;
	free(obj->norm_x);
	free(obj->norm_f);
	free(obj->trans_x);
	free(obj->basic_f);
	free(obj->temp_x1);
	free(obj->temp_x2);
	free(obj->temp_x3);
	free(obj->temp_x4);
	free(obj->weight);
	free(obj->sigma);
	free(obj->lam);
	free(obj->bias);
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		free(obj->l[i][j]);
	}
	for (i = 0; i < obj->nfunc; i++) {
		free(obj->o[i]);
		free(obj->l[i]);
	}
	for (i = 0; i < obj->nreal; i++) {
		free(obj->g[i]);
	}
	free(obj->o);
	free(obj->l);
	free(obj->g);
	if (obj->nfunc == 5) {
		// TODO: free Af5
		// TODO: free Bf5
	} else if (obj->nfunc == 12) {
		// TODO: free Af12
		// TODO: free Bf12
		// TODO: free alphaf12
	}
	free(obj);
	obj = NULL;
	return;
}

static char *dirname;

void set_directory(char *dir) {
	dirname = dir;
}

FILE *myopen(const char *file, const char *mode) {
	FILE *fpt;
	fpt = fopen(file, mode);
	if (fpt == NULL) {
		char fname[500];
		sprintf(fname, "%s/%s", dirname, file);
		fpt = fopen(fname, mode);
	}
	return fpt;
}

void initialize_f1(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("sphere_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -450.0;
	return;
}

void initialize_f2(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("schwefel_102_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -450.0;
	return;
}

void initialize_f3(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	if (obj->nreal == 2)
		fpt = myopen("elliptic_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("elliptic_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("elliptic_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("elliptic_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("high_cond_elliptic_rot_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -450.0;
	return;
}

void initialize_f4(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("schwefel_102_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -450.0;
	return;
}

void initialize_f5(CEC2005data *obj) {
	int i, j, index;
	char c;
	FILE *fpt;
	obj->Af5 = (long double **)malloc(obj->nreal * sizeof(long double));
	for (i = 0; i < obj->nreal; i++) {
		obj->Af5[i] = (long double *)malloc(obj->nreal * sizeof(long double));
	}
	obj->Bf5 = (long double *)malloc(obj->nreal * sizeof(long double));
	fpt = myopen("schwefel_206_data.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	for (i = 0; i < obj->nreal; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->Af5[i][j]));
			//printf("\n A[%d][%d] = %LE",i+1,j+1,A[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	if (obj->nreal % 4 == 0) {
		index = obj->nreal / 4;
	} else {
		index = obj->nreal / 4 + 1;
	}
	for (i = 0; i < index; i++) {
		obj->o[0][i] = -100.0;
	}
	index = (3 * obj->nreal) / 4 - 1;
	for (i = index; i < obj->nreal; i++) {
		obj->o[0][i] = 100.0;
	}
	for (i = 0; i < obj->nreal; i++) {
		obj->Bf5[i] = 0.0;
		for (j = 0; j < obj->nreal; j++) {
			obj->Bf5[i] += obj->Af5[i][j] * obj->o[0][j];
		}
	}
	obj->bias[0] = -310.0;
	return;
}

void initialize_f6(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("rosenbrock_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		obj->o[i][j] -= 1.0;
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = 390.0;
	return;
}

void initialize_f7(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	if (nreal == 2)
		fpt = myopen("griewank_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("griewank_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("griewank_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("griewank_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("griewank_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	bias[0] = -180.0;
	return;
}

void initialize_f8(CEC2005data *obj) {
	int i, j;
	int index;
	FILE *fpt;
	if (obj->nreal == 2)
		fpt = myopen("ackley_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("ackley_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("ackley_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("ackley_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("ackley_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	index = obj->nreal / 2;
	for (i = 1; i <= index; i++) {
		obj->o[0][2 * i - 2] = -32.0;
	}
	obj->bias[0] = -140.0;
	return;
}

void initialize_f9(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("rastrigin_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -330.0;
	return;
}

void initialize_f10(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	if (obj->nreal == 2)
		fpt = myopen("rastrigin_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("rastrigin_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("rastrigin_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("rastrigin_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nreal; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("rastrigin_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &o[i][j]);
		//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -330.0;
	return;
}

void initialize_f11(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	if (obj->nreal == 2)
		fpt = myopen("weierstrass_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("weierstrass_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("weierstrass_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("weierstrass_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("weierstrass_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = 90.0;
	return;
}

void initialize_f12(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	char c;
	obj->Af12 = (long double **)malloc(obj->nreal * sizeof(long double));
	obj->Bf12 = (long double **)malloc(obj->nreal * sizeof(long double));
	obj->alphaf12 = (long double *)malloc(obj->nreal * sizeof(long double));
	for (i = 0; i < obj->nreal; i++) {
		Af12[i] = (long double *)malloc(obj->nreal * sizeof(long double));
		Bf12[i] = (long double *)malloc(obj->nreal * sizeof(long double));
	}
	fpt = myopen("schwefel_213_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) {
		for (j = 0; j < obj->nreal; j++)	{
			fscanf(fpt, "%Lf", &(obj->Af12[i][j]));
			//            printf("\n Af12[%d][%d] = %LE",i+1,j+1,Af12[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	if (i != 100) {
		for (i = obj->nreal; i < 100; i++) {
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	for (i = 0; i < obj->nreal; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->Bf12[i][j]));
			//            printf("\n B[%d][%d] = %LE",i+1,j+1,B[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	if (i != 100) {
		for (i = obj->nreal; i < 100; i++) {
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	/* Reading alpha */
	for (i = 0; i < obj->nreal; i++) {
		fscanf(fpt, "%Lf", &(obj->alphaf12[i]));
		//        printf("\n alpha[%d] = %LE",i+1,alpha[i]);
	}
	fclose(fpt);
	obj->bias[0] = -460.0;
	return;
}

void initialize_f13(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	fpt = myopen("EF8F2_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			obj->o[i][j] -= 1.0;
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
	}
	fclose(fpt);
	obj->bias[0] = -130.0;
	return;
}

void initialize_f14(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	if (obj->nreal == 2)
		fpt = myopen("E_ScafferF6_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("E_ScafferF6_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("E_ScafferF6_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("E_ScafferF6_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nreal; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->g[i][j]));
		//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("E_ScafferF6_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &(obj->o[i][j]));
		//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	obj->bias[0] = -300.0;
	return;
}

void initialize_f15(CEC2005data *obj) {
	int i, j;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	obj->lam[0] = 1.0;
	obj->lam[1] = 1.0;
	obj->lam[2] = 10.0;
	obj->lam[3] = 10.0;
	obj->lam[4] = 1.0 / 12.0;
	obj->lam[5] = 1.0 / 12.0;
	obj->lam[6] = 5.0 / 32.0;
	obj->lam[7] = 5.0 / 32.0;
	obj->lam[8] = 1.0 / 20.0;
	obj->lam[9] = 1.0 / 20.0;
	obj->global_bias = 120.0;
	return;
}

void initialize_f16(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func1_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func1_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func1_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func1_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	obj->lam[0] = 1.0;
	obj->lam[1] = 1.0;
	obj->lam[2] = 10.0;
	obj->lam[3] = 10.0;
	obj->lam[4] = 1.0 / 12.0;
	obj->lam[5] = 1.0 / 12.0;
	obj->lam[6] = 5.0 / 32.0;
	obj->lam[7] = 5.0 / 32.0;
	obj->lam[8] = 1.0 / 20.0;
	obj->lam[9] = 1.0 / 20.0;
	obj->global_bias = 120.0;
	return;
}

void initialize_f17(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func1_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func1_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func1_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func1_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	obj->lam[0] = 1.0;
	obj->lam[1] = 1.0;
	obj->lam[2] = 10.0;
	obj->lam[3] = 10.0;
	obj->lam[4] = 1.0 / 12.0;
	obj->lam[5] = 1.0 / 12.0;
	obj->lam[6] = 5.0 / 32.0;
	obj->lam[7] = 5.0 / 32.0;
	obj->lam[8] = 1.0 / 20.0;
	obj->lam[9] = 1.0 / 20.0;
	obj->global_bias = 120.0;
	return;
}

void initialize_f18(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			// printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	for (j = 0; j < obj->nreal; j++) {
		obj->o[9][j] = 0.0;
	}
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	for (i = 0; i < obj->nreal; i++) {
		obj->o[obj->nfunc - 1][i] = 0.0;
	}
	obj->sigma[0] = 1.0;
	obj->sigma[1] = 2.0;
	obj->sigma[2] = 1.5;
	obj->sigma[3] = 1.5;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 1.0;
	obj->sigma[6] = 1.5;
	obj->sigma[7] = 1.5;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 5.0 / 16.0;
	obj->lam[1] = 5.0 / 32.0;
	obj->lam[2] = 2.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 1.0 / 10.0;
	obj->lam[5] = 1.0 / 20.0;
	obj->lam[6] = 20.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 6.0;
	obj->lam[9] = 1.0 / 12.0;
	obj->global_bias = 10.0;
	return;
}

void initialize_f19(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	for (j = 0; j < obj->nreal; j++) {
		obj->o[9][j] = 0.0;
	}
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	for (i = 0; i < obj->nreal; i++) {
		obj->o[obj->nfunc - 1][i] = 0.0;
	}
	obj->sigma[0] = 0.1;
	obj->sigma[1] = 2.0;
	obj->sigma[2] = 1.5;
	obj->sigma[3] = 1.5;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 1.0;
	obj->sigma[6] = 1.5;
	obj->sigma[7] = 1.5;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 0.5 / 32.0;
	obj->lam[1] = 5.0 / 32.0;
	obj->lam[2] = 2.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 1.0 / 10.0;
	obj->lam[5] = 1.0 / 20.0;
	obj->lam[6] = 20.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 6.0;
	obj->lam[9] = 1.0 / 12.0;
	obj->global_bias = 10.0;
	return;
}

void initialize_f20(CEC2005data *obj) {
	int i, j, k, index;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	/**
	 * Daniel: Arreglado errror
	 */
	for (j = 0; j < obj->nreal; j++) {
		obj->o[9][j] = 0.0;
	}
	index = obj->nreal / 2;
	for (i = 1; i <= index; i++) {
		obj->o[0][2 * i - 1] = 5.0;
	}
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	for (i = 0; i < obj->nreal; i++) {
		obj->o[nfunc - 1][i] = 0.0;
	}
	obj->sigma[0] = 1.0;
	obj->sigma[1] = 2.0;
	obj->sigma[2] = 1.5;
	obj->sigma[3] = 1.5;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 1.0;
	obj->sigma[6] = 1.5;
	obj->sigma[7] = 1.5;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 5.0 / 16.0;
	obj->lam[1] = 5.0 / 32.0;
	obj->lam[2] = 2.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 1.0 / 10.0;
	obj->lam[5] = 1.0 / 20.0;
	obj->lam[6] = 20.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 6.0;
	obj->lam[9] = 1.0 / 12.0;
	obj->global_bias = 10.0;
	return;
}

void initialize_f21(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func3_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func3_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func3_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func3_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	obj->sigma[0] = 1.0;
	obj->sigma[1] = 1.0;
	obj->sigma[2] = 1.0;
	obj->sigma[3] = 1.0;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 2.0;
	obj->sigma[6] = 2.0;
	obj->sigma[7] = 2.0;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 1.0 / 4.0;
	obj->lam[1] = 1.0 / 20.0;
	obj->lam[2] = 5.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 5.0;
	obj->lam[5] = 1.0;
	obj->lam[6] = 50.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 8.0;
	obj->lam[9] = 1.0 / 40.0;
	obj->global_bias = 360.0;
	return;
}

void initialize_f22(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func3_HM_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func3_HM_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func3_HM_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func3_HM_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	obj->sigma[0] = 1.0;
	obj->sigma[1] = 1.0;
	obj->sigma[2] = 1.0;
	obj->sigma[3] = 1.0;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 2.0;
	obj->sigma[6] = 2.0;
	obj->sigma[7] = 2.0;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 1.0 / 4.0;
	obj->lam[1] = 1.0 / 20.0;
	obj->lam[2] = 5.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 5.0;
	obj->lam[5] = 1.0;
	obj->lam[6] = 50.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 8.0;
	obj->lam[9] = 1.0 / 40.0;
	obj->global_bias = 360.0;
	return;
}

void initialize_f23(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func3_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func3_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func3_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func3_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	obj->sigma[0] = 1.0;
	obj->sigma[1] = 1.0;
	obj->sigma[2] = 1.0;
	obj->sigma[3] = 1.0;
	obj->sigma[4] = 1.0;
	obj->sigma[5] = 2.0;
	obj->sigma[6] = 2.0;
	obj->sigma[7] = 2.0;
	obj->sigma[8] = 2.0;
	obj->sigma[9] = 2.0;
	obj->lam[0] = 1.0 / 4.0;
	obj->lam[1] = 1.0 / 20.0;
	obj->lam[2] = 5.0;
	obj->lam[3] = 1.0;
	obj->lam[4] = 5.0;
	obj->lam[5] = 1.0;
	obj->lam[6] = 50.0;
	obj->lam[7] = 10.0;
	obj->lam[8] = 1.0 / 8.0;
	obj->lam[9] = 1.0 / 40.0;
	obj->global_bias = 360.0;
	return;
}

void initialize_f24(CEC2005data *obj) {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func4_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			fscanf(fpt, "%Lf", &(obj->o[i][j]));
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	fclose(fpt);
	if (obj->nreal == 2)
		fpt = myopen("hybrid_func4_M_D2.txt", "r");
	if (obj->nreal == 10)
		fpt = myopen("hybrid_func4_M_D10.txt", "r");
	if (obj->nreal == 30)
		fpt = myopen("hybrid_func4_M_D30.txt", "r");
	if (obj->nreal == 50)
		fpt = myopen("hybrid_func4_M_D50.txt", "r");
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			for (k = 0; k < obj->nreal; k++) {
				fscanf(fpt, "%Lf", &(obj->l[i][j][k]));
				//printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do {
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	for (i = 0; i < 10; i++) {
		obj->sigma[i] = 2.0;
	}
	obj->lam[0] = 10.0;
	obj->lam[1] = 1.0 / 4.0;
	obj->lam[2] = 1.0;
	obj->lam[3] = 5.0 / 32.0;
	obj->lam[4] = 1.0;
	obj->lam[5] = 1.0 / 20.0;
	obj->lam[6] = 1.0 / 10.0;
	obj->lam[7] = 1.0;
	obj->lam[8] = 1.0 / 20.0;
	obj->lam[9] = 1.0 / 20.0;
	obj->global_bias = 260.0;
	return;
}

void initialize_f25(CEC2005data *obj) {
	initialize_f24(obj);
}

void initialize(CEC2005data *obj) {
	int num = obj->nfunc;
	if (num == 1) {
		initialize_f1(obj);
	} else if (num == 2) {
		initialize_f2(obj);
	} else if (num == 3) {
		initialize_f3(obj);
	} else if (num == 4) {
		initialize_f4(obj);
	} else if (num == 5) {
		initialize_f5(obj);
	} else if (num == 6) {
		initialize_f6(obj);
	} else if (num == 7) {
		initialize_f7(obj);
	} else if (num == 8) {
		initialize_f8(obj);
	} else if (num == 9) {
		initialize_f9(obj);
	} else if (num == 10) {
		initialize_f10(obj);
	} else if (num == 11) {
		initialize_f11(obj);
	} else if (num == 12) {
		initialize_f12(obj);
	} else if (num == 13) {
		initialize_f13(obj);
	} else if (num == 14) {
		initialize_f14(obj);
	} else if (num == 15) {
		initialize_f15(obj);
		calc_benchmark_norm_f15(obj);
	} else if (num == 16) {
		initialize_f16(obj);
		calc_benchmark_norm_f16(obj);
	} else if (num == 17) {
		initialize_f17(obj);
		calc_benchmark_norm_f17(obj);
	} else if (num == 18) {
		initialize_f18(obj);
		calc_benchmark_norm_f18(obj);
	} else if (num == 19) {
		initialize_f19(obj);
		calc_benchmark_norm_f19(obj);
	} else if (num == 20) {
		initialize_f20(obj);
		calc_benchmark_norm_f20(obj);
	} else if (num == 21) {
		initialize_f21(obj);
		calc_benchmark_norm_f21(obj);
	} else if (num == 22) {
		initialize_f22(obj);
		calc_benchmark_norm_f22(obj);
	} else if (num == 23) {
		initialize_f23(obj);
		calc_benchmark_norm_f23(obj);
	} else if (num == 24) {
		initialize_f24(obj);
		calc_benchmark_norm_f24(obj);
	} else if (num == 25) {
		initialize_f25(obj);
		calc_benchmark_norm_f25(obj);
	} else {
		printf("Error: num %d no v�lido\n", num);
		exit(1);
	}
}

long double calc_benchmark_func(long double *x, CEC2005data *fdata) {
	int num = fdata->nfunc;
	if (num == 1) {
		return calc_benchmark_func_f1(x, fdata);
	} else if (num == 2) {
		return calc_benchmark_func_f2(x, fdata);
	} else if (num == 3) {
		return calc_benchmark_func_f3(x, fdata);
	} else if (num == 4) {
		return calc_benchmark_func_f4(x, fdata);
	} else if (num == 5) {
		return calc_benchmark_func_f5(x, fdata);
	} else if (num == 6) {
		return calc_benchmark_func_f6(x, fdata);
	} else if (num == 7) {
		return calc_benchmark_func_f7(x, fdata);
	} else if (num == 8) {
		return calc_benchmark_func_f8(x, fdata);
	} else if (num == 9) {
		return calc_benchmark_func_f9(x, fdata);
	} else if (num == 10) {
		return calc_benchmark_func_f10(x, fdata);
	} else if (num == 11) {
		return calc_benchmark_func_f11(x, fdata);
	} else if (num == 12) {
		return calc_benchmark_func_f12(x, fdata);
	} else if (num == 13) {
		return calc_benchmark_func_f13(x, fdata);
	} else if (num == 14) {
		return calc_benchmark_func_f14(x, fdata);
	} else if (num == 15) {
		return calc_benchmark_func_f15(x, fdata);
	} else if (num == 16) {
		return calc_benchmark_func_f16(x, fdata);
	} else if (num == 17) {
		return calc_benchmark_func_f17(x, fdata);
	} else if (num == 18) {
		return calc_benchmark_func_f18(x, fdata);
	} else if (num == 19) {
		return calc_benchmark_func_f19(x, fdata);
	} else if (num == 20) {
		return calc_benchmark_func_f20(x, fdata);
	} else if (num == 21) {
		return calc_benchmark_func_f21(x, fdata);
	} else if (num == 22) {
		return calc_benchmark_func_f22(x, fdata);
	} else if (num == 23) {
		return calc_benchmark_func_f23(x, fdata);
	} else if (num == 24) {
		return calc_benchmark_func_f24(x, fdata);
	} else if (num == 25) {
		return calc_benchmark_func_f25(x, fdata);
	} else {
		printf("Error: num %d no v�lido\n", num);
		exit(1);
	}
	return 0;
}

bool isBound_cec2005(void) {
	if (nfunc != 7 && nfunc != 25)     {
		return true;
	} else {
		return false;
	}
}

/**
 * fun: Number/index of benchmark function
 * dim: Number of variable/components of the selected benchmark functions
 */
CEC2005data* init_cec2005(int fun, int dim) {
	// Asigno valor a las variables
	int nfunc = fun, nreal = dim;
	if (nfunc < 1) {
		fprintf(stderr, "\n Wrong value of 'nfun' entered\n");
		exit(0);
	}
	if (nreal != 2 && nreal != 10 && nreal != 30 && nreal != 50) {
		fprintf(stderr, "\n Wrong value of 'nreal' entered, only 2, 10, 30, 50 variables are supported\n");
		exit(0);
	}
	CEC2005data *fdata = allocate_memory(nreal, nfunc);
	initialize(fdata);
	return fdata;
}

void finish_cec2005(CEC2005data *fdata) {
	free_memory(fdata);
}

struct cec2005function {
	int ident;       /**< Par�metro identificador de la funci�n */
	char* name;     /**< Nombre de la funci�n (para poder mostrarla en pantalla */
	double range[2]; /**< Guarda los valores m�nimos y m�ximos para dicha funci�n */
	double optime;   /**< Valor �ptimo */
};

typedef struct cec2005function CECFUNCTION;

static CECFUNCTION cec2005Fun[] = {
	{1, "Shifted Sphere Function", {-100, 100}, -450},
	{2, "Shifted Schwefel's Problem 1.2", {-100, 100}, -450},
	{3, "Shifted Rotated High Conditioned Elliptic Function", {-100, 100}, -450},
	{4, "Shifted Schwefel's Problem 1.2 with Noise in Fitness", {-100, 100}, -450},
	{5, "Schwefel's Problem 2.6 with Global Optimum on Bounds", {-100, 100}, -310},
	{6, "Shifted Rosenbrock's Function", {-100, 100}, 390},
	{7, "Shifted Rotated Griewnk's Function without Bounds", {0, 600}, -180},
	{8, "Shifted Rotated Ackley's Function with Global Optimum on Bounds", {-32, 32}, -140},
	{9, "Shifted Rastrigin's Function", {-5, 5}, -330},
	{10, "Shifted Rotated Rastrigin's Function", {-5, 5}, -330},
	{11, "Shifted Rotated Weierstrass Function", {-0.5, 0.5}, 90},
	{12, "Schwefel's Problem 2.13", {-PI, PI}, -460},
	{13, "Expanded Functions", {-5, 5}, -130},
	{14, "Shifted Rotated Expanded Scaffer's F6 Function", {-100, 100}, -300},
	{15, "Hybrid Composition Function 1", {-5, 5}, 120},
	{16, "Rotated Hybrid Composition Function 1", {-5, 5}, 120},
	{17, "Rotated Hybrid Composition Function 1 with Noise in Fitness", {-5, 5}, 120},
	{18, "Rotated Hybrid Composition Function 2", {-5, 5}, 10},
	{19, "Rotated Hybrid Composition Function 2 with Narrow Basin Global Optimum", {-5, 5}, 10},
	{20, "Rotated Hybrid Composition Function 2 with Global Optimum on the Bounds", {-5, 5}, 10},
	{21, "Rotated Hybrid Composition Function 3", {-5, 5}, 360},
	{22, "Rotated Hybrid Composition Function 3 with High Condition Number Matrix", {-5, 5}, 360},
	{23, "Non-Continuous Rotated Hybrid Composition Function 3", {-5, 5}, 360},
	{24, "Rotated Hybrid Composition Function 4", {-5, 5}, 260},
	{25, "Rotated Hybrid Composition Function 4 without Bounds", {2, 5}, 260}
};

static int cec2005FunSize = 25;

/**
 * Permite obtener la informaci�n sobre la funci�n
 */
void getInfo_cec2005(int fun, char *name, double *min, double *max, double *optime) {
	int id;
	assert(fun >= 0 && fun <= cec2005FunSize);
	id = fun - 1;
	strcpy(name, cec2005Fun[id].name);
	*min = cec2005Fun[id].range[0];
	*max = cec2005Fun[id].range[1];
	*optime = cec2005Fun[id].optime;
}

/**
 * La funci�n eval del CEC2005
 */
double eval_cec2005_ld(const long double *x, int ndim) {
	double optime = cec2005Fun[nfunc - 1].optime;
	long double *y = (long double *)x;
	double fit = calc_benchmark_func(y) - optime;
	assert(fit >= 0);
	return fit;
}

double eval_cec2005(const double *x, int ndim) {
	double optime = cec2005Fun[nfunc - 1].optime;
	long double y[ndim];
	for (int i = 0; i < ndim; ++i) {
		y[i] = x[i];
	}
	double fit = calc_benchmark_func(y) - optime;
	if (fit < 0) {
		fprintf(stderr, "Value: %le\tOptime: %le\n", fit, optime);
	}
	assert(fit >= 0);
	return fit;
}

int main() {
	return 0;
}
