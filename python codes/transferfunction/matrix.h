#pragma once
#ifndef MATRIX__H
#define MATRIX__H

#include "mbed.h"
#include <vector>
#include <cmath>

using std::vector;



class Matrix {

public:
	//    Matrix();
	Matrix(unsigned int, unsigned int);
	double& operator()(const unsigned int &, const unsigned int &);
	void print(void);
	unsigned row_size(void);
	unsigned column_size(void);
	Matrix operator+(double);
	Matrix operator-(double);
	Matrix operator*(double);
	Matrix operator/(double);
	Matrix operator+(Matrix &);
	Matrix operator-(Matrix &);
	Matrix operator*(Matrix &);
	Matrix trans(void);
	Matrix find_minor(unsigned int, unsigned int);
	double det(void);
	Matrix inv(void);

private:
	unsigned _row_size;
	unsigned _column_size;
	vector<vector<double> > _matrix;

};

#endif /* MATRIX__H) */