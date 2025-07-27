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
	norm_x = (long double *)malloc(nreal * sizeof(long double));
	norm_f = (long double *)malloc(nfunc * sizeof(long double));
	trans_x = (long double *)malloc(nreal * sizeof(long double));
	basic_f = (long double *)malloc(nfunc * sizeof(long double));
	temp_x1 = (long double *)malloc(nreal * sizeof(long double));
	temp_x2 = (long double *)malloc(nreal * sizeof(long double));
	temp_x3 = (long double *)malloc(nreal * sizeof(long double));
	temp_x4 = (long double *)malloc(nreal * sizeof(long double));
	weight = (long double *)malloc(nfunc * sizeof(long double));
	sigma = (long double *)malloc(nfunc * sizeof(long double));
	lambda = (long double *)malloc(nfunc * sizeof(long double));
	bias = (long double *)malloc(nfunc * sizeof(long double));
	o = (long double **)malloc(nfunc * sizeof(long double));
	l = (long double ***)malloc(nfunc * sizeof(long double));
	g = (long double **)malloc(nreal * sizeof(long double));
	for (i = 0; i < nfunc; i++) {
		o[i] = (long double *)malloc(nreal * sizeof(long double));
		l[i] = (long double **)malloc(nreal * sizeof(long double));
	}
	for (i = 0; i < nreal; i++) {
		g[i] = (long double *)malloc(nreal * sizeof(long double));
	}
	for (i = 0; i < nfunc; i++) for (j = 0; j < nreal; j++) {
		l[i][j] = (long double *)malloc(nreal * sizeof(long double));

	}
	/* Do some trivial (common) initialization here itself */
	C = 2000.0;
	for (i = 0; i < nreal; i++) {
		norm_x[i] = 5.0;
		trans_x[i] = 0.0;
		temp_x1[i] = 0.0;
		temp_x2[i] = 0.0;
		temp_x3[i] = 0.0;
		temp_x4[i] = 0.0;
		for (j = 0; j < nreal; j++) {
			if (i == j) {
				g[i][j] = 1.0;
			} else {
				g[i][j] = 0.0;
			}
		}
	}
	for (i = 0; i < nfunc; i++) {
		basic_f[i] = 0.0;
		norm_f[i] = 0.0;
		weight[i] = 1.0 / (long double)nfunc;
		sigma[i] = 1.0;
		lambda[i] = 1.0;
		bias[i] = 100.0 * (long double)i;
		for (j = 0; j < nreal; j++) {
			o[i][j] = 0.0;
			for (k = 0; k < nreal; k++) {
				if (j == k) {
					l[i][j][k] = 1.0;
				} else {
					l[i][j][k] = 0.0;
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
	free(obj->lambda);
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
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			// FIXME
			fscanf(fpt, "%Lf", &o[i][j]);
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
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
	for (i = 0; i < obj->nfunc; i++) {
		for (j = 0; j < obj->nreal; j++) {
			// FIXME
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
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
		// FIXME
		fscanf(fpt, "%Lf", &g[i][j]);
		//printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("high_cond_elliptic_rot_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < obj->nfunc; i++) for (j = 0; j < obj->nreal; j++) {
		fscanf(fpt, "%Lf", &o[i][j]);
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
		// FIXME
		fscanf(fpt, "%Lf", &o[i][j]);
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
			// FIXME
			fscanf(fpt, "%Lf", &o[i][j]);
			//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do {
			// FIXME
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	for (i = 0; i < obj->nreal; i++) {
		for (j = 0; j < obj->nreal; j++) {
			// FIXME
			fscanf(fpt, "%Lf", &Af5[i][j]);
			//printf("\n A[%d][%d] = %LE",i+1,j+1,A[i][j]);
		}
		do {
			// FIXME
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
		fscanf(fpt, "%Lf", &o[i][j]);
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
	for (i = 0; i < nreal; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &g[i][j]);
		//printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("griewank_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &o[i][j]);
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
	if (nreal == 2)
		fpt = myopen("ackley_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("ackley_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("ackley_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("ackley_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nreal; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &g[i][j]);
		//printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
	}
	fclose(fpt);
	fpt = myopen("ackley_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &o[i][j]);
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	index = nreal / 2;
	for (i = 1; i <= index; i++) {
		o[0][2 * i - 2] = -32.0;
	}
	bias[0] = -140.0;
	return;
}

void initialize_f9() {
	int i, j;
	FILE *fpt;
	fpt = myopen("rastrigin_func_data.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++) for (j = 0; j < nreal; j++) {
		fscanf(fpt, "%Lf", &o[i][j]);
		//printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
	}
	fclose(fpt);
	bias[0] = -330.0;
	return;
}

void initialize_f10() {
	int i, j;
	FILE *fpt;
	if (nreal == 2)
		fpt = myopen("rastrigin_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("rastrigin_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("rastrigin_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("rastrigin_M_D50.txt", "r");
	if (fpt == NULL) {
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nreal; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &g[i][j]);
			//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
		}
	}
	fclose(fpt);
	fpt = myopen("rastrigin_func_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
	}
	fclose(fpt);
	bias[0] = -330.0;
	return;
}

void initialize_f11() {
	int i, j;
	FILE *fpt;
	if (nreal == 2)
		fpt = myopen("weierstrass_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("weierstrass_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("weierstrass_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("weierstrass_M_D50.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nreal; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &g[i][j]);
			//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
		}
	}
	fclose(fpt);
	fpt = myopen("weierstrass_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
	}
	fclose(fpt);
	bias[0] = 90.0;
	return;
}

void initialize_f12() {
	int i, j;
	FILE *fpt;
	char c;
	Af12 = (long double **)malloc(nreal * sizeof(long double));
	Bf12 = (long double **)malloc(nreal * sizeof(long double));
	alphaf12 = (long double *)malloc(nreal * sizeof(long double));
	for (i = 0; i < nreal; i++)
	{
		Af12[i] = (long double *)malloc(nreal * sizeof(long double));
		Bf12[i] = (long double *)malloc(nreal * sizeof(long double));
	}
	fpt = myopen("schwefel_213_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	/* Reading A */
	for (i = 0; i < nreal; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &Af12[i][j]);
			//            printf("\n Af12[%d][%d] = %LE",i+1,j+1,Af12[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	if (i != 100)
	{
		for (i = nreal; i < 100; i++)
		{
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	/* Reading B */
	for (i = 0; i < nreal; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &Bf12[i][j]);
			//            printf("\n B[%d][%d] = %LE",i+1,j+1,B[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
	}
	if (i != 100)
	{
		for (i = nreal; i < 100; i++)
		{
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
	}
	/* Reading alpha */
	for (i = 0; i < nreal; i++)
	{
		fscanf(fpt, "%Lf", &alphaf12[i]);
		//        printf("\n alpha[%d] = %LE",i+1,alpha[i]);
	}
	fclose(fpt);
	bias[0] = -460.0;
	return;
}

void initialize_f13() {
	int i, j;
	FILE *fpt;
	fpt = myopen("EF8F2_func_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			o[i][j] -= 1.0;
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
	}
	fclose(fpt);
	bias[0] = -130.0;
	return;
}

void initialize_f14() {
	int i, j;
	FILE *fpt;
	if (nreal == 2)
		fpt = myopen("E_ScafferF6_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("E_ScafferF6_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("E_ScafferF6_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("E_ScafferF6_M_D50.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nreal; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &g[i][j]);
			//            printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]);
		}
	}
	fclose(fpt);
	fpt = myopen("E_ScafferF6_func_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
	}
	fclose(fpt);
	bias[0] = -300.0;
	return;
}

void initialize_f15() {
	int i, j;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	lambda[0] = 1.0;
	lambda[1] = 1.0;
	lambda[2] = 10.0;
	lambda[3] = 10.0;
	lambda[4] = 1.0 / 12.0;
	lambda[5] = 1.0 / 12.0;
	lambda[6] = 5.0 / 32.0;
	lambda[7] = 5.0 / 32.0;
	lambda[8] = 1.0 / 20.0;
	lambda[9] = 1.0 / 20.0;
	global_bias = 120.0;
	return;
}

void initialize_f16() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func1_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func1_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func1_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func1_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	lambda[0] = 1.0;
	lambda[1] = 1.0;
	lambda[2] = 10.0;
	lambda[3] = 10.0;
	lambda[4] = 1.0 / 12.0;
	lambda[5] = 1.0 / 12.0;
	lambda[6] = 5.0 / 32.0;
	lambda[7] = 5.0 / 32.0;
	lambda[8] = 1.0 / 20.0;
	lambda[9] = 1.0 / 20.0;
	global_bias = 120.0;
	return;
}

void initialize_f17() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func1_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func1_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func1_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func1_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func1_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	lambda[0] = 1.0;
	lambda[1] = 1.0;
	lambda[2] = 10.0;
	lambda[3] = 10.0;
	lambda[4] = 1.0 / 12.0;
	lambda[5] = 1.0 / 12.0;
	lambda[6] = 5.0 / 32.0;
	lambda[7] = 5.0 / 32.0;
	lambda[8] = 1.0 / 20.0;
	lambda[9] = 1.0 / 20.0;
	global_bias = 120.0;
	return;
}

void initialize_f18() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);

	/**
	 * Daniel: Arreglado errror
	 */
	for (j = 0; j < nreal; j++)
	{
		o[9][j] = 0.0;
	}

	if (nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	for (i = 0; i < nreal; i++)
	{
		o[nfunc - 1][i] = 0.0;
	}
	sigma[0] = 1.0;
	sigma[1] = 2.0;
	sigma[2] = 1.5;
	sigma[3] = 1.5;
	sigma[4] = 1.0;
	sigma[5] = 1.0;
	sigma[6] = 1.5;
	sigma[7] = 1.5;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 5.0 / 16.0;
	lambda[1] = 5.0 / 32.0;
	lambda[2] = 2.0;
	lambda[3] = 1.0;
	lambda[4] = 1.0 / 10.0;
	lambda[5] = 1.0 / 20.0;
	lambda[6] = 20.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 6.0;
	lambda[9] = 1.0 / 12.0;
	global_bias = 10.0;
	return;
}

void initialize_f19() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);

	/**
	 * Daniel: Arreglado errror
	 */
	for (j = 0; j < nreal; j++)
	{
		o[9][j] = 0.0;
	}

	if (nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	for (i = 0; i < nreal; i++)
	{
		o[nfunc - 1][i] = 0.0;
	}
	sigma[0] = 0.1;
	sigma[1] = 2.0;
	sigma[2] = 1.5;
	sigma[3] = 1.5;
	sigma[4] = 1.0;
	sigma[5] = 1.0;
	sigma[6] = 1.5;
	sigma[7] = 1.5;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 0.5 / 32.0;
	lambda[1] = 5.0 / 32.0;
	lambda[2] = 2.0;
	lambda[3] = 1.0;
	lambda[4] = 1.0 / 10.0;
	lambda[5] = 1.0 / 20.0;
	lambda[6] = 20.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 6.0;
	lambda[9] = 1.0 / 12.0;
	global_bias = 10.0;
	return;
}

void initialize_f20() {
	int i, j, k;
	int index;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func2_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);

	/**
	 * Daniel: Arreglado errror
	 */
	for (j = 0; j < nreal; j++)
	{
		o[9][j] = 0.0;
	}

	index = nreal / 2;
	for (i = 1; i <= index; i++)
	{
		o[0][2 * i - 1] = 5.0;
	}
	if (nreal == 2)
		fpt = myopen("hybrid_func2_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func2_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func2_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func2_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	for (i = 0; i < nreal; i++)
	{
		o[nfunc - 1][i] = 0.0;
	}
	sigma[0] = 1.0;
	sigma[1] = 2.0;
	sigma[2] = 1.5;
	sigma[3] = 1.5;
	sigma[4] = 1.0;
	sigma[5] = 1.0;
	sigma[6] = 1.5;
	sigma[7] = 1.5;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 5.0 / 16.0;
	lambda[1] = 5.0 / 32.0;
	lambda[2] = 2.0;
	lambda[3] = 1.0;
	lambda[4] = 1.0 / 10.0;
	lambda[5] = 1.0 / 20.0;
	lambda[6] = 20.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 6.0;
	lambda[9] = 1.0 / 12.0;
	global_bias = 10.0;
	return;
}

void initialize_f21() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func3_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func3_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func3_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func3_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	sigma[0] = 1.0;
	sigma[1] = 1.0;
	sigma[2] = 1.0;
	sigma[3] = 1.0;
	sigma[4] = 1.0;
	sigma[5] = 2.0;
	sigma[6] = 2.0;
	sigma[7] = 2.0;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 1.0 / 4.0;
	lambda[1] = 1.0 / 20.0;
	lambda[2] = 5.0;
	lambda[3] = 1.0;
	lambda[4] = 5.0;
	lambda[5] = 1.0;
	lambda[6] = 50.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 8.0;
	lambda[9] = 1.0 / 40.0;
	global_bias = 360.0;
	return;
}

void initialize_f22() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func3_HM_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func3_HM_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func3_HM_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func3_HM_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//  printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	sigma[0] = 1.0;
	sigma[1] = 1.0;
	sigma[2] = 1.0;
	sigma[3] = 1.0;
	sigma[4] = 1.0;
	sigma[5] = 2.0;
	sigma[6] = 2.0;
	sigma[7] = 2.0;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 1.0 / 4.0;
	lambda[1] = 1.0 / 20.0;
	lambda[2] = 5.0;
	lambda[3] = 1.0;
	lambda[4] = 5.0;
	lambda[5] = 1.0;
	lambda[6] = 50.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 8.0;
	lambda[9] = 1.0 / 40.0;
	global_bias = 360.0;
	return;
}

void initialize_f23() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func3_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func3_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func3_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func3_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func3_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
				/*printf("\n got here \n");*/
			} while (c != '\n');
		}
		//        printf("\n");
	}
	sigma[0] = 1.0;
	sigma[1] = 1.0;
	sigma[2] = 1.0;
	sigma[3] = 1.0;
	sigma[4] = 1.0;
	sigma[5] = 2.0;
	sigma[6] = 2.0;
	sigma[7] = 2.0;
	sigma[8] = 2.0;
	sigma[9] = 2.0;
	lambda[0] = 1.0 / 4.0;
	lambda[1] = 1.0 / 20.0;
	lambda[2] = 5.0;
	lambda[3] = 1.0;
	lambda[4] = 5.0;
	lambda[5] = 1.0;
	lambda[6] = 50.0;
	lambda[7] = 10.0;
	lambda[8] = 1.0 / 8.0;
	lambda[9] = 1.0 / 40.0;
	global_bias = 360.0;
	return;
}

void initialize_f24() {
	int i, j, k;
	FILE *fpt;
	char c;
	fpt = myopen("hybrid_func4_data.txt", "r");
	if (fpt == NULL)
	{
		fprintf(stderr, "\n Error: Cannot open input file for reading \n");
		exit(0);
	}
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			fscanf(fpt, "%Lf", &o[i][j]);
			//            printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]);
		}
		do
		{
			fscanf(fpt, "%c", &c);
		} while (c != '\n');
		//        printf("\n");
	}
	fclose(fpt);
	if (nreal == 2)
		fpt = myopen("hybrid_func4_M_D2.txt", "r");
	if (nreal == 10)
		fpt = myopen("hybrid_func4_M_D10.txt", "r");
	if (nreal == 30)
		fpt = myopen("hybrid_func4_M_D30.txt", "r");
	if (nreal == 50)
		fpt = myopen("hybrid_func4_M_D50.txt", "r");
	for (i = 0; i < nfunc; i++)
	{
		for (j = 0; j < nreal; j++)
		{
			for (k = 0; k < nreal; k++)
			{
				fscanf(fpt, "%Lf", &l[i][j][k]);
				//                printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]);
			}
			do
			{
				fscanf(fpt, "%c", &c);
			} while (c != '\n');
		}
		//        printf("\n");
	}
	for (i = 0; i < 10; i++)
	{
		sigma[i] = 2.0;
	}
	lambda[0] = 10.0;
	lambda[1] = 1.0 / 4.0;
	lambda[2] = 1.0;
	lambda[3] = 5.0 / 32.0;
	lambda[4] = 1.0;
	lambda[5] = 1.0 / 20.0;
	lambda[6] = 1.0 / 10.0;
	lambda[7] = 1.0;
	lambda[8] = 1.0 / 20.0;
	lambda[9] = 1.0 / 20.0;
	global_bias = 260.0;
	return;
}

void initialize_f25() {
	initialize_f24();
}

void initialize(void) {
	int num = nfunc;
	if (num == 1) {
		initialize_f1();
	} else if (num == 2) {
		initialize_f2();
	} else if (num == 3) {
		initialize_f3();
	} else if (num == 4) {
		initialize_f4();
	} else if (num == 5) {
		initialize_f5();
	} else if (num == 6) {
		initialize_f6();
	} else if (num == 7) {
		initialize_f7();
	} else if (num == 8) {
		initialize_f8();
	} else if (num == 9) {
		initialize_f9();
	} else if (num == 10) {
		initialize_f10();
	} else if (num == 11) {
		initialize_f11();
	} else if (num == 12) {
		initialize_f12();
	} else if (num == 13) {
		initialize_f13();
	} else if (num == 14) {
		initialize_f14();
	} else if (num == 15) {
		initialize_f15();
		calc_benchmark_norm_f15();
	} else if (num == 16) {
		initialize_f16();
		calc_benchmark_norm_f16();
	} else if (num == 17) {
		initialize_f17();
		calc_benchmark_norm_f17();
	} else if (num == 18) {
		initialize_f18();
		calc_benchmark_norm_f18();
	} else if (num == 19) {
		initialize_f19();
		calc_benchmark_norm_f19();
	} else if (num == 20) {
		initialize_f20();
		calc_benchmark_norm_f20();
	} else if (num == 21) {
		initialize_f21();
		calc_benchmark_norm_f21();
	} else if (num == 22) {
		initialize_f22();
		calc_benchmark_norm_f22();
	} else if (num == 23) {
		initialize_f23();
		calc_benchmark_norm_f23();
	} else if (num == 24) {
		initialize_f24();
		calc_benchmark_norm_f24();
	} else if (num == 25) {
		initialize_f25();
		calc_benchmark_norm_f25();
	} else {
		printf("Error: num %d no v�lido\n", num);
		exit(1);
	}
}

typedef long double (*tBenchmark)(long double *);

static tBenchmark bench; 

tBenchmark set_calc_benchmark_func(int num) {
	if (num == 1) {
		return calc_benchmark_func_f1;
	}
	else if (num == 2) {
		return calc_benchmark_func_f2;
	}
	else if (num == 3) {
		return calc_benchmark_func_f3;
	}
	else if (num == 4) {
		return calc_benchmark_func_f4;
	}
	else if (num == 5) {
		return calc_benchmark_func_f5;
	}
	else if (num == 6) {
		return calc_benchmark_func_f6;
	}
	else if (num == 7) {
		return calc_benchmark_func_f7;
	}
	else if (num == 8) {
		return calc_benchmark_func_f8;
	}
	else if (num == 9) {
		return calc_benchmark_func_f9;
	}
	else if (num == 10) {
		return calc_benchmark_func_f10;
	}
	else if (num == 11) {
		return calc_benchmark_func_f11;
	}
	else if (num == 12) {
		return calc_benchmark_func_f12;
	}
	else if (num == 13) {
		return calc_benchmark_func_f13;
	}
	else if (num == 14) {
		return calc_benchmark_func_f14;
	}
	else if (num == 15) {
		return calc_benchmark_func_f15;
	}
	else if (num == 16) {
		return calc_benchmark_func_f16;
	}
	else if (num == 17) {
		return calc_benchmark_func_f17;
	}
	else if (num == 18) {
		return calc_benchmark_func_f18;
	}
	else if (num == 19) {
		return calc_benchmark_func_f19;
	}
	else if (num == 20) {
		return calc_benchmark_func_f20;
	}
	else if (num == 21) {
		return calc_benchmark_func_f21;
	}
	else if (num == 22) {
		return calc_benchmark_func_f22;
	}
	else if (num == 23) {
		return calc_benchmark_func_f23;
	}
	else if (num == 24) {
		return calc_benchmark_func_f24;
	}
	else if (num == 25) {
		return calc_benchmark_func_f25;
	}
	else {
		printf("Error: num %d no v�lido\n", num);
		exit(1);
	}

	return NULL;
}

long double calc_benchmark_func(long double *x) {
	static int notinit = 1;

	if (notinit) {
		bench = set_calc_benchmark_func(nfunc);
	}

	return (*bench)(x);
}

bool isBound_cec2005(void) {
	if (nfunc != 7 && nfunc != 25)     {
		return true;
	} else {
		return false;
	}
}

void init_cec2005(int fun, int dim) {
	// Asigno valor a las variables
	nfunc = fun;
	nreal = dim;

	if (nfunc < 1) {
		fprintf(stderr, "\n Wrong value of 'nfun' entered\n");
		exit(0);
	}
	if (nreal != 2 && nreal != 10 && nreal != 30 && nreal != 50) {
		fprintf(stderr, "\n Wrong value of 'nreal' entered, only 2, 10, 30, 50 variables are supported\n");
		exit(0);
	}

	/* nreal and nfun need to be initialized before calling these routines */
	/* Routine to allocate memory to global variables */
	allocate_memory();

	/* Routine the initalize global variables */
	initialize();
}

void finish_cec2005(void) {
	free_memory();
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
