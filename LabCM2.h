#pragma once
#include<iostream>
#include<math.h>

using namespace std;


double kx(double x);

double qx(double x);

double fx(double x);

double IntegralSimpsKX(double a, double b, int N);

double IntegralSimpsQX(double a, double b, int N);

double IntegralSimpsFX(double a, double b, int N);

double Jx(double h, double x, double x1, double x2, double n);

double Dx(double h, double x, double x1, double x2, double n);

double Ax(double h, double x, double x1, double x2, double n);

double* progonka(int n, double* a, double* b, double* c, double* f);

double* Difference(int n, double* x);

void Osnov();