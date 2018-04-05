#ifndef MATHE_H
#define MATHE_H

#define PI 3.14159265359
#define mean 0.0
#define Dr 5e-3    //rotational diffusion rate
#include<cstdbool>

extern double sigma;
extern double N[10];
extern int    *cluster;
extern bool   check_cal_cluster;

double random_number(int iseed);
double Gaussian_noise(double x1, double x2, int iseed);// x1 = mean  x2 = standard deviation
double array_mean(int *a, int len);        // calculate the mean of array
void   cal_std_deviation(void);
void   cal_cluster_size(void);
void   cal_recursive(int p1, int n, bool *visited, int **nl);

#endif
