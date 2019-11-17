#include<iostream>
#include "LabCm2.h"


using namespace std;


double kx(double x)
{
	if ((x >= 0) && (x <= 1))
	{
		if ((x >= 0) && (x < 0.5))
			return pow((x + 1), 2);
		if ((x >= 0.5) && x <= 1)
			return (pow(x, 2));
	}
	else
		cout << " sm kx" << endl;
}

double qx(double x)
{
	double e = 2.718281828459;
	if ((x >= 0) && (x <= 1))
	{
		if ((x >= 0) && (x < 0.5))
			return pow(e, (-1) * x) * sqrt(e);
		if ((x >= 0.5) && x <= 1)
			return (pow(e, x) * (0.6065306597));
	}
	else 
		cout << " sm qx" << endl;
}

double fx(double x)
{
	double e = 2.718281828459;
	//cout << x << endl;
	if ((x >= 0) && (x <= 1))
	{
		if ((x >= 0) && (x < 0.5))
			return cos((3.141592653589)*x);
		if ((x >= 0.5) && x <= 1)
			return (pow(e, x) * (0.6065306597));
	}
	else 
		cout << " sm fx" << endl;
}

double IntegralSimpsKX(double a, double b, int N)
{
	double S, xx, h;
	S = 0.0;
	h = (b - a) / (double)N;
	xx = a + h;
	while (xx < b)
	{
		S += 4.0 * (1.0/kx(xx));
		xx += h;
		if (xx >= b)
			break;
		S += 2.0 * (1.0 / kx(xx));
		xx += h;
	}
	S += (1.0 / kx(a)) + (1.0 / kx(b));
	S *= h / 3.0;
	return S;
}

double IntegralSimpsQX(double a, double b, int N)
{
	double S , xx, h;
	S = 0.0;
	h = (b - a) /(double) N;
	xx = a + h;
	while (xx < b)
	{
		S += 4.0 * qx(xx);
		xx += h;
		if (xx >= b)
			break;
		S += 2.0 * qx(xx);
		xx += h;
	}
	S += qx(a) + qx(b);
	S *= h / 3.0;
	return S;
}

double IntegralSimpsFX(double a, double b, int N)
{
	double S , xx, h;
	S = 0.0;
	h = (b - a) / (double)N;
	xx = a + h;
	while (xx < b)
	{
		S += 4.0 * fx(xx);
		xx += h;
		if (xx >= b)
			break;
		S += 2.0 * fx(xx);
		xx += h;
	}
	S += fx(a) + fx(b);
	S *= h / 3.0;
	return S;
}

double Jx(double h, double x, double x1, double x2, int n)
{
	if ((x2 <= 0.5) || (x1 >= 0.5))
		return fx(x);
	else
		return ((1.0/h)*(IntegralSimpsFX(x1,0.5, n) + IntegralSimpsFX(0.5, x2, n)));
}

double Dx(double h, double x, double x1, double x2, int n)
{
	if ((x2 <= 0.5) || (x1 >= 0.5))
		return qx(x);
	else
		return ((1.0 / h) * (IntegralSimpsQX(x1,0.5,n) + IntegralSimpsQX(0.5,x2,n)));
}

double Ax(double h, double x, double x1, double x2, int n)
{
	if ((x2 <= 0.5) || (x1 >= 0.5))
	{
		return pow(kx(x), 2);
	}
		
	else
	{
		return (h)*(1.0/(IntegralSimpsKX(x1,0.5,n) + IntegralSimpsKX(0.5, x2,  n)));
	}
}

double* progonka(int n, double* a, double* b, double* c, double* f)
{
	//n - число уравнений (строк матрицы)
	//b - диагональ, лежащая над главной 
	//c - главная диагональ матрицы 
	//a - диагональ, лежащая под главной 
	//f - правая часть 
	//y - решение

	double *m1, *m2;
	m1  = new double[n];
	m2  = new double[n];
	m1[0] = 0; m2[0] = 0;

	double* y;
	y = new double[10000];
	for (int i = 1; i < n; i++)
	{
		m1[i] = b[i - 1] / (c[i - 1] - m1[i - 1]*a[i-1]);
		m2[i] = (m1[i - 1] + a[i - 1] * b[i - 1]) / (c[i - 1] - m1[i - 1] * a[i - 1]);
	}

	y[n - 1] = (f[n - 1] );

	for (int i = n - 1; i >= 0; i--)
		y[i] = (f[i] - b[i] * y[i + 1]) / c[i];

	return y;
}

double* Difference(int n, double *x)
{
	double h = 1.0 / (double)n;
	//Формирование матриц для метода прогонки

	double* a; // Диагональ под главной диагональю
	a = new double[10000];
	a[0] = 0;
	for (int i = 1; i < n; i++)
	{
		a[i] = 1.0 ;
	}
	

	double* b;// Диагональ над главной диагональю
	b = new double[10000];
	b[n - 1] = 0;
	for (int i = 0; i < n - 1; i++)
	{
		double j = (double)i;
		b[i] = (-1.0)*Ax(h, x[i], (j) * h, (j + 1.0)* h, n);
	}

	double* c; // Главная диагональ
	c = new double[10000];
	c[0] = 1;
	c[n - 1] = 1;
	for (int i = 1; i < n - 1; i++)
	{
		c[i] = Ax(h, x[i], (i)* h, (i + 1.0) * h, n) + pow(h,2)* Dx(h, x[i], (i - 0.5) * h, (i + 0.5) * h, n) + 1.0;
	}

	double* f;
	f = new double[10000];
	for (int i = 0; i < n - 1; i++)
	{
		f[i] = pow(h,2)*Jx(h, x[i + 1], (i + 0.5) * h, (i + 1.5) * h, n);
	}

	double *rezult;
	rezult = new double[10000];
	rezult = progonka(n, a, b, c, f);

	return rezult;
}

void Osnov()
{
	int n;

	cout << "Enrer n" << endl;
	cin >> n;

	double h; 
	h = 1 /(double) n;

	double *x, *x12;
	x = new double[n];
	x12 = new double[n];
	for (int i = 0; i < n ; i++) // узлы сеток 
	{
		double xi = i * h;
		x[i] = xi;
		//cout << xi << endl;
	}
	for (int i = 0; i < n ; i++)
	{
		x12[i] = (i + 0.5) * h;
		//cout << x12[i] << endl;
	}

	double *v;
	double* v2;
	v = new double[n];
	v2 = new double[n];
	v = Difference(n, x);
	v2 = Difference(n, x12);

	double *e2;
	e2 = new double[n];
	int b = 0;
	for (int i = 0; i < n; i++)
	{
		e2[i] = v2[i] - v[i];
		if (e2[i] < 0)
			e2[i] = (-1) * e2[i];
		
		if (e2[i] < 0.0000005)
			b = 1;
		else b = 0;
	}

	for (int i = 0; i < n; i++)
	{
		cout << "Number node" << i << endl;
		cout << "x:  " << x[i] << "          v(x):    " << v[i] << "       v2(x):    " << v2[i] << endl;
		cout << "v(x) - v2(x):   " << e2[i];
		if (b == 1) cout << "      True" << endl;
		else cout <<"        False" << endl;
		cout << endl;
		cout << endl;
	}
}

