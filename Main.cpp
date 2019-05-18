#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "Matrix.h"
#include <ctime>
#include <stdexcept>
#include <numeric>

using namespace std;

const double a = 1;
const double c1 = 1.5;
const double c2 = -3;

const long double tMin = 0.0;
const long double tMax = 1.0;
const long double xMin = 0.0;
const long double xMax = 1.0;

const int xSteps = 5;
const int tSteps = 500;
const long double xStep = 1.0 / xSteps - 1;
const long double tStep = 1.0 / tSteps -1 ;

unsigned int exact_time;
unsigned int approx_time;

double function_w(double previousW, double currentW, double nextW, double q, double h) {
	return currentW + q*(pow(currentW, -1.0)*((nextW - 2 * currentW + previousW) / pow(h, 2.0)) - pow(currentW, -2.0)*pow((nextW - previousW) / 2 * h, 2.0));
}

double exactFunction(double x, double t) {
	auto val = (long double)1 / (c1*x - pow(c1, 2.0)*t + c2);
	if (fabs(val) >= 10000000000)
		throw std::logic_error("watafak");
	return val;
}

void assessErrors(int xSteps, int tSteps, Matrix exact, Matrix approx);
void saveToFiles(string methodName, Matrix exact, Matrix approx, int xSteps, int tSteps);
void solve();

void solve() {
	Matrix exact = Matrix(tSteps + 1, xSteps + 1);
	Matrix approx = Matrix(tSteps + 1, xSteps + 1);
	int i = 0; // for t
	int j = 0; // for x
	
	long double t;
	long double x;

	unsigned int exact_start_time = clock();
	for (i = 0, t = tMin, x = xMin; i <= tSteps; ++i, t += tStep) {
		for (j = 0; j <= xSteps; ++j, x += xStep) {
			exact(i, j) = exactFunction(x, t);
		}
		x = xMin;
	}
	unsigned int exact_end_time = clock();
	unsigned int exact_time = exact_end_time - exact_start_time;


	unsigned int approx_start_time = clock();
	
	for (i = 0, t = tMin; i <= tSteps; ++i, t += tStep) {
		approx(i, 0) = exactFunction(0, t);
		approx(i, xSteps) = exactFunction(xSteps, t);
	}
	
	for (j = 0, x = xMin; j <= xSteps; ++j, x += xStep) {
		approx(0, j) = exactFunction(x, t);
	}
	
	for (i = 1; i <= tSteps; ++i) {
#pragma omp parallel for
		for (j = 1; j < xSteps; ++j) {
			approx(i, j) = function_w(approx(i - 1, j - 1), approx(i - 1, j), approx(i - 1, j + 1), tStep, xStep);
		}
	}
	unsigned int approx_end_time = clock();
	unsigned int approx_time = approx_end_time - approx_start_time;

	assessErrors(xSteps, tSteps, exact, approx);
	saveToFiles("Matrix", exact, approx, xSteps, tSteps);
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
			long double relativeError = absoluteError / exact(i, j) * 100.0;
			totalAbsoluteError += absoluteError;
			totalRelativeError += relativeError;
			maxAbsoluteError = max(maxAbsoluteError, absoluteError);
			maxRelativeError = max(maxRelativeError, relativeError);
		}
	}
	long double averageAbsoluteError = totalAbsoluteError / (tSteps * xSteps);
	long double averageRelativeError = totalRelativeError / (tSteps * xSteps);


	cout << "Середня абсолютна похибка : " << averageAbsoluteError << "\n"
	     << "Середня вiдносна похибка : " << averageRelativeError << "\n"
		 << "Mаксимальна абсолютна похибка : " << maxAbsoluteError << "\n"
		 << "Максимальна вiдносна похибка : " << maxRelativeError << "\n\n"
	     << "Час виконання точного розв'язку : " << exact_time << "\n"
		 << "Час виконання знайденого розв'язку : " << approx_time << endl;

	for (int i = 0; i <= tSteps; ++i) {
		for (int j = 0; j <= xSteps; ++j) {
			cout << fixed << setprecision(5) << approx(i, j) << "\t";
		}
		cout << endl;
	}

	cout << endl;
	
	for (int i = 0; i <= tSteps; ++i) {
		for (int j = 0; j <= xSteps; ++j) {
			cout << fixed << setprecision(5) << exact(i, j) << "\t";
		}
		cout << endl;
	}

}

void saveToFiles(string methodName, Matrix exact, Matrix approx, int xSteps, int tSteps) {
	ofstream fout1;
	ofstream fout2("approx.txt", ios_base::out);
	fout1.open("exact.txt", ios_base::out);
	for (int i = 0; i <= tSteps; ++i) {
		for (int j = 0; j <= xSteps; ++j) {
			fout1 << fixed << setprecision(5) << exact(i, j) << ",";
			fout2 << fixed << setprecision(5) << approx(i, j) << ",";
		}
		fout1 << endl;
		fout2 << endl;
	}
	fout1.close();
	fout2.close();
}

int main(void) {
	setlocale(LC_ALL, "Russian");
	omp_set_num_threads(4);
	solve();
	system("pause");
	return 0;
}