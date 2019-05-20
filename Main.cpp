#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "Matrix.h"
#include <ctime>


using namespace std;
int val;

const double a = 1;
const double c1 = 0.3;
const double c2 = 1;

const long double tMin = 0.0;
const long double tMax = 1.0;
const long double xMin = 0.0;
const long double xMax = 1.0;

const int xSteps = 8;
const int tSteps = 1000
;
const long double xStep = 1.0 / (xSteps -1) ;
const long double tStep = 1.0 / (tSteps - 1);

double function_w(double previousW, double currentW, double nextW, double q, double h) {
	return currentW + a*q*(((1/currentW)*((nextW - 2 * currentW + previousW) / pow(h, 2.0))) - (pow(currentW, -2.0)*pow((nextW - previousW) / 2 * h, 2.0)));
}

double exactFunction(double x, double t) {
	return (long double)1 / (c1*x - pow(c1, 2.0)*t + c2);
}

void assessErrors(int xSteps, int tSteps, Matrix exact, Matrix approx);
void saveToFiles(Matrix exact, Matrix approx, int xSteps, int tSteps);
void solve();

void solve() {
	Matrix exact = Matrix(tSteps + 1, xSteps + 1);
	Matrix approx = Matrix(tSteps + 1, xSteps + 1);
	int i = 0; // for t
	int j = 0; // for x
	
	long double t;
	long double x;

	cout << "Виконати послiдовно(0) або паралельно(1)? ";
	cin >> val;
	unsigned int start_time_exact = clock();
    #pragma omp parallel if (val)
	if (omp_in_parallel()) {
		for (i = 0, t = tMin, x = xMin; i <= tSteps; ++i, t += tStep) {
			#pragma omp parallel for
			for (j = 0; j <= xSteps; ++j, x += xStep) {
				exact(i, j) = exactFunction(x, t);
			}
			x = xMin;
		}
		cout << "Час виконання точного розв'язку паралельно: " << clock() - start_time_exact << endl;
	}
	else {
		for (i = 0, t = tMin, x = xMin; i <= tSteps; ++i, t += tStep) {
			for (j = 0; j <= xSteps; ++j, x += xStep) {
				exact(i, j) = exactFunction(x, t);
			}
			x = xMin;
		}
		cout << "Час виконання точного розв'язку послiдовно: " << clock() - start_time_exact << endl;
	}
	
	


	unsigned int start_time_approx = clock();
	
	for (i = 0, t = tMin; i <= tSteps; ++i, t += tStep) {
		approx(i, 0) = exactFunction(0, t);
		approx(i, xSteps) = exactFunction(xSteps, t);
	}
	
	for (j = 0, x = xMin; j <= xSteps; ++j, x += xStep) {
		approx(0, j) = exactFunction(x, t);
	}

#pragma omp parallel if (val)
	if (omp_in_parallel()) {
		for (i = 1; i <= tSteps; ++i) {
		#pragma omp parallel for
			for (j = 1; j < xSteps; ++j) {
				approx(i, j) = function_w(approx(i - 1, j - 1), approx(i - 1, j), approx(i - 1, j + 1), tStep, xStep);
			}
		}
		cout << "Час виконання наближеного розв'язку паралельно: " << clock() - start_time_approx << endl;
	}
	else {
		for (i = 1; i <= tSteps; ++i) {
			for (j = 1; j < xSteps; ++j) {
				approx(i, j) = function_w(approx(i - 1, j - 1), approx(i - 1, j), approx(i - 1, j + 1), tStep, xStep);
			}
		}
		cout << "Час виконання наближеного розв'язку послiдовно: " << clock() - start_time_approx << endl;
	}
	

	assessErrors(xSteps, tSteps, exact, approx);
	saveToFiles(exact, approx, xSteps, tSteps);
	return ;
}

void assessErrors(int xSteps, int tSteps, Matrix exact, Matrix approx) {
	long double totalAbsoluteError = 0.0;
	long double totalRelativeError = 0.0;
	long double maxAbsoluteError = 0.0;
	long double maxRelativeError = 0.0;
	for (int i = 1; i <= tSteps; ++i) {
		for (int j = 1; j < xSteps; ++j) {
			long double absoluteError = abs(approx(i, j) - exact(i, j));
			long double relativeError = abs(absoluteError / exact(i, j) * 100.0);
			totalAbsoluteError += absoluteError;
			totalRelativeError += relativeError;
			maxAbsoluteError = max(maxAbsoluteError, absoluteError);
			maxRelativeError = max(maxRelativeError, relativeError);
		}
	}
	long double averageAbsoluteError = totalAbsoluteError / (tSteps * xSteps);
	long double averageRelativeError = totalRelativeError / (tSteps * xSteps);


	cout << "Середня абсолютна похибка : " << fixed << setprecision(5) << averageAbsoluteError << endl;
	cout << "Середня вiдносна похибка : " << fixed << setprecision(5) << averageRelativeError << "%" << endl;
	cout << "Mаксимальна абсолютна похибка : " << fixed << setprecision(5) << maxAbsoluteError << endl;
	cout << "Максимальна вiдносна похибка : " << fixed << setprecision(5) << maxRelativeError << "%\n" << endl;

}

void saveToFiles(Matrix exact, Matrix approx, int xSteps, int tSteps) {
	ofstream fout1("exact.txt");
	ofstream fout2("approx.txt");
	
	for (int i = 0; i <= tSteps; ++i) {
		for (int j = 0; j <= xSteps; ++j) {
			fout1 << fixed << setprecision(5) << exact(i, j) << "\t";
			fout2 << fixed << setprecision(5) << approx(i, j) << "\t";
		}
		fout1 << endl;
		fout2 << endl;
	}
	fout1.close();
	fout2.close();
}

int main(void) {
	setlocale(LC_ALL, "Russian");
	solve();
	system("pause");
	return 0;
}