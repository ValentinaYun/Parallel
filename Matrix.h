#pragma once
class Matrix {
private:
	long double *arr;
	int rows;
	int cols;
public:
	Matrix(int rows, int cols);
	long double &operator()(const int &row, const int &col);
	long double *getAsArray();
};