#include "benchfunctions.h"

#include <math.h>

double calc_benchmark_func_f1(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_sphere (fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f2(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_schwefel (fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f3(double *x, CEC2005data *fdata) {
	int i;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	double y = 0.0;
	for (i = 0; i < fdata->nreal; i++) {
		y += fdata->trans_x[i] * fdata->trans_x[i] * pow(1.0e6, i / (fdata->nreal - 1.0));
	}
	return y + fdata->bias[0];
}

double calc_benchmark_func_f4(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	double y = calc_schwefel(fdata->trans_x, fdata->nreal) * (1.0 + 0.4 * fabs(randomnormaldeviate()));
	return y + fdata->bias[0];
}

double calc_benchmark_func_f5(double *x, CEC2005data *fdata) {
	int i, j;
    double y;
	fdata->basic_f[0] = -INF;
	for (i =0 ; i < fdata->nreal; i++) {
		y = 0.0;
		for (j = 0; j < fdata->nreal; j++) {
			y += fdata->Af5[i][j] * x[j];
		}
		y = fabs(y - fdata->Bf5[i]);
		if (fdata->basic_f[0] < y) {
			fdata->basic_f[0] = y;
		}
	}
	return fdata->basic_f[0] + fdata->bias[0];
}

double calc_benchmark_func_f6(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_rosenbrock(fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f7(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_griewank(fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f8(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_ackley(fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f9(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_rastrigin(fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f10(double *x, CEC2005data *fdata) {
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	return calc_rastrigin(fdata->trans_x, fdata->nreal) + fdata->bias[0];
}

double calc_benchmark_func_f11(double *x, CEC2005data *fdata) {
	int i;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	double y = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	return y + fdata->bias[0];
}

double calc_benchmark_func_f12(double *x, CEC2005data *fdata) {
	int i, j;
	double sum1, sum2, y = 0.0;
	for (i = 0; i < fdata->nreal; i++) {
		sum1 = 0.0, sum2 = 0.0;
		for (j = 0; j < fdata->nreal; j++) {
			sum1 += fdata->Af12[i][j] * sin(fdata->alphaf12[j]) + fdata->Bf12[i][j] * cos(fdata->alphaf12[j]);
			sum2 += fdata->Af12[i][j] * sin(x[j]) + fdata->Bf12[i][j] * cos(x[j]);
		}
		y += pow((sum1 - sum2), 2.0);
	}
	return y + fdata->bias[0];
}

double calc_benchmark_func_f13(double *x, CEC2005data *fdata) {
	int i;
	double tmp, y = 0.0;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	for (i = 0; i < fdata->nreal - 1; i++) {
		tmp = 100.0 * pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		y += (tmp * tmp) / 4000.0 - cos(tmp) + 1.0;
	}
	tmp = 100.0 * pow((fdata->trans_x[fdata->nreal - 1] * fdata->trans_x[fdata->nreal - 1] - fdata->trans_x[0]), 2.0) + 1.0 * pow((fdata->trans_x[fdata->nreal - 1] - 1.0), 2.0);
	y += (tmp * tmp) / 4000.0 - cos(tmp) + 1.0 + fdata->bias[0];
	return y;
}

double calc_benchmark_func_f14(double *x, CEC2005data *fdata) {
	int i;
	double temp1, temp2, y = 0.0;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	for (i = 0; i < fdata->nreal-1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		y += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	y += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0)) + fdata->bias[0];
	return y;
}

void calc_benchmark_norm_f15(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f15(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double y;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i  < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	y = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		y += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return y;
}

void calc_benchmark_norm_f16(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f16(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double y;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	y = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		y += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return y;
}

void calc_benchmark_norm_f17(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f17(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double y = 0.0;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	for (i = 0; i < nfunc; i++) {
		y += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	y *= (1.0 + 0.2 * fabs(randomnormaldeviate()));
	y += fdata->global_bias;
	return y;
}

void calc_benchmark_norm_f18(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f18(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double res;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f19(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f19(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double res;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i = 0; i < nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f20(CEC2005data *fdata) {
	int i;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f20(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double res;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_sphere(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_sphere(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i=0; i<nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f21(CEC2005data *fdata) {
	int i;
	double temp1, temp2, temp;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->norm_f[0] += 0.5 + (temp1 - 0.5) / pow(temp2, 2.0);
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[fdata->nreal-1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->norm_f[0] += 0.5 + (temp1 - 0.5) / pow(temp2, 2.0);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->norm_f[1] += 0.5 + (temp1 - 0.5) / pow(temp2, 2.0);
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->norm_f[1] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = 0.0;
	for (i = 0; i < fdata->nreal-1; i++) {
		temp = 100.0 * pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->norm_f[4] += (temp * temp)/ 4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal - 1] * fdata->trans_x[fdata->nreal - 1] - fdata->trans_x[0]), 2.0) + 1.0 * pow((fdata->trans_x[fdata->nreal - 1] - 1.0), 2.0);
	fdata->norm_f[4] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0 * pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->norm_f[5] += (temp * temp) /4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0 * pow((fdata->trans_x[fdata->nreal - 1] * fdata->trans_x[fdata->nreal - 1] - fdata->trans_x[0]), 2.0) + 1.0 * pow((fdata->trans_x[fdata->nreal - 1] - 1.0), 2.0);
	fdata->norm_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f21(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double temp1, temp2, temp, res;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[0] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->basic_f[0] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[1] += 0.5 + (temp1 - 0.5) / (pow(temp2, 2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0 * pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->basic_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f22(CEC2005data *fdata) {
	int i;
	double temp1, temp2, temp;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0)))),2.0);
		temp2 = 1.0 + 0.001*(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0));
		fdata->norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0)))),2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0));
	fdata->norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0)))),2.0);
		temp2 = 1.0 + 0.001*(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0));
		fdata->norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0)))),2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0));
	fdata->norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal - 1] * fdata->trans_x[fdata->nreal - 1] - fdata->trans_x[0]), 2.0) + 1.0 * pow((fdata->trans_x[fdata->nreal - 1] - 1.0), 2.0);
	fdata->norm_f[5] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f22(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double temp1, temp2, temp, res;
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = 0.0;
	for (i = 0; i< fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1],2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = 0.0;
	for (i = 0; i < fdata->nreal-1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1],2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0)))),2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0));
	fdata->basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0 * pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f23(CEC2005data *fdata) {
	int i;
	double temp1, temp2, temp;
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0)))),2.0);
		temp2 = 1.0 + 0.001*(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0));
		fdata->norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0)))),2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0));
	fdata->norm_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0)))),2.0);
		temp2 = 1.0 + 0.001*(pow(fdata->trans_x[i],2.0)+pow(fdata->trans_x[i+1],2.0));
		fdata->norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0)))),2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal-1],2.0)+pow(fdata->trans_x[0],2.0));
	fdata->norm_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->norm_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->norm_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	return;
}

double calc_benchmark_func_f23(double *x, CEC2005data *fdata) {
	int i, a, nfunc = 10;
	double temp1, temp2, temp, b, res;
	for (i = 0; i < fdata->nreal; i++) {
		if (fabs(x[i] - fdata->o[0][i]) >= 0.5) {
			res = 2.0 * x[i];
			a = (int) res;
			b = fabs(res - a);
			if (b < 0.5) {
				fdata->temp_x4[i] = a / 2.0;
			} else {
				if (res <= 0.0) {
					fdata->temp_x4[i] = (a - 1.0) / 2.0;
				} else {
					fdata->temp_x4[i] = (a + 1.0) / 2.0;
				}
			}
		} else {
			fdata->temp_x4[i] = x[i];
		}
	}
	transform (fdata->temp_x4, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->basic_f[0] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform (fdata->temp_x4, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp1 = pow((sin(sqrt(pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1],2.0)))), 2.0);
		temp2 = 1.0 + 0.001 * (pow(fdata->trans_x[i], 2.0) + pow(fdata->trans_x[i + 1], 2.0));
		fdata->basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	}
	temp1 = pow((sin(sqrt(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0)))), 2.0);
	temp2 = 1.0 + 0.001*(pow(fdata->trans_x[fdata->nreal - 1], 2.0) + pow(fdata->trans_x[0], 2.0));
	fdata->basic_f[1] += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
	transform (fdata->temp_x4, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (fdata->temp_x4, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (fdata->temp_x4, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0 * pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal - 1] * fdata->trans_x[fdata->nreal - 1] - fdata->trans_x[0]), 2.0) + 1.0 * pow((fdata->trans_x[fdata->nreal - 1]-1.0), 2.0);
	fdata->basic_f[4] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform (fdata->temp_x4, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i]*fdata->trans_x[i]-fdata->trans_x[i+1]),2.0) + 1.0*pow((fdata->trans_x[i]-1.0),2.0);
		fdata->basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[5] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (fdata->temp_x4, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (fdata->temp_x4, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (fdata->temp_x4, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (fdata->temp_x4, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = calc_griewank(fdata->trans_x, fdata->nreal);
	for (i = 0; i < fdata->nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(fdata->temp_x4, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f24(CEC2005data *fdata) {
	int i;
	double temp;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform_norm (0, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[0] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform_norm (1, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		fdata->norm_f[1] += ExpandedF6(fdata->trans_x[i], fdata->trans_x[i + 1]);
	}
	fdata->norm_f[1] += ExpandedF6(fdata->trans_x[fdata->nreal - 1], fdata->trans_x[0]);
	transform_norm (2, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[2] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i]-1.0),2.0);
		fdata->norm_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->norm_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform_norm (3, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[3] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform_norm (4, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[4] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform_norm (5, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform_norm (6, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[6] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		fdata->norm_f[6] += nc_schaffer(fdata->trans_x[i], fdata->trans_x[i + 1]);
	}
	fdata->norm_f[6] += nc_schaffer(fdata->trans_x[fdata->nreal - 1], fdata->trans_x[0]);
	transform_norm (7, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[7] = nc_rastrigin(fdata->trans_x, fdata->nreal, fdata->temp_x4);
	transform_norm (8, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[8] = 0.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_f[8] += fdata->trans_x[i] * fdata->trans_x[i] * pow(1.0e6, i / (fdata->nreal - 1.0));
	}
	transform_norm (9, fdata->nreal, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->g, fdata->l);
	fdata->norm_f[9] = calc_sphere(fdata->trans_x, fdata->nreal) * (1.0 + 0.1 * fabs(randomnormaldeviate()));
	return;
}

double calc_benchmark_func_f24(double *x, CEC2005data *fdata) {
	int i, nfunc = 10;
	double temp, res;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->norm_x[i] = 0.0;
	}
	transform (x, 0, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[0] = calc_weierstrass(fdata->trans_x, fdata->nreal) - calc_weierstrass(fdata->norm_x, fdata->nreal);
	transform (x, 1, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[1] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		fdata->basic_f[1] += ExpandedF6(fdata->trans_x[i], fdata->trans_x[i + 1]);
	}
	fdata->basic_f[1] += ExpandedF6(fdata->trans_x[fdata->nreal - 1], fdata->trans_x[0]);
	transform (x, 2, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[2] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		temp = 100.0*pow((fdata->trans_x[i] * fdata->trans_x[i] - fdata->trans_x[i + 1]), 2.0) + 1.0 * pow((fdata->trans_x[i] - 1.0), 2.0);
		fdata->basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	temp = 100.0*pow((fdata->trans_x[fdata->nreal-1]*fdata->trans_x[fdata->nreal-1]-fdata->trans_x[0]),2.0) + 1.0*pow((fdata->trans_x[fdata->nreal-1]-1.0),2.0);
	fdata->basic_f[2] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	transform (x, 3, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[3] = calc_ackley(fdata->trans_x, fdata->nreal);
	transform (x, 4, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[4] = calc_rastrigin(fdata->trans_x, fdata->nreal);
	transform (x, 5, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[5] = calc_griewank(fdata->trans_x, fdata->nreal);
	transform (x, 6, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[6] = 0.0;
	for (i = 0; i < fdata->nreal - 1; i++) {
		fdata->basic_f[6] += nc_schaffer(fdata->trans_x[i], fdata->trans_x[i+1]);
	}
	fdata->basic_f[6] += nc_schaffer(fdata->trans_x[fdata->nreal - 1], fdata->trans_x[0]);
	transform (x, 7, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[7] = nc_rastrigin(fdata->trans_x, fdata->nreal, fdata->temp_x4);
	transform (x, 8, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[8] = 0.0;
	for (i = 0; i < fdata->nreal; i++) {
		fdata->basic_f[8] += fdata->trans_x[i] * fdata->trans_x[i] * pow(1.0e6, i / (fdata->nreal - 1.0));
	}
	transform (x, 9, fdata->nreal, fdata->temp_x1, fdata->temp_x2, fdata->temp_x3, fdata->trans_x, fdata->lam, fdata->o, fdata->g, fdata->l);
	fdata->basic_f[9] = (calc_sphere(fdata->trans_x, fdata->nreal)) * (1.0 + 0.1* fabs(randomnormaldeviate()));
	for (i = 0; i < nfunc; i++) {
		fdata->basic_f[i] *= fdata->C / fdata->norm_f[i];
	}
	calc_weight(x, nfunc, fdata->nreal, fdata->weight, fdata->sigma, fdata->o);
	res = fdata->global_bias;
	for (i = 0; i < nfunc; i++) {
		res += fdata->weight[i] * (fdata->basic_f[i] + fdata->bias[i]);
	}
	return (res);
}

void calc_benchmark_norm_f25(CEC2005data *fdata) {
	calc_benchmark_norm_f24(fdata);
}

double calc_benchmark_func_f25(double *x, CEC2005data *fdata) {
	return calc_benchmark_func_f24(x, fdata);   
}
