#include <random>
#include <stdio.h>
#include "Matrix.h"
Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols) {
	arr = new long double[rows * cols];
}
long double &Matrix::operator()(const int &row, const int &col) {
	return this -> arr[row * cols + col];
}

long double *Matrix::getAsArray() {
	return this -> arr;
}