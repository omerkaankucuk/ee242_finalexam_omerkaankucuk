#include "matrix.h"

Matrix::Matrix(unsigned int row_size, unsigned int column_size) {
	_row_size = row_size;
	_column_size = column_size;
	_matrix.resize(_row_size);
	for (int i = 0; i < _row_size; i++) {
		_matrix[i].resize(_column_size);
	}
}

void Matrix::print(void) {
	printf("Matrix:\n");
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			printf("%g     ", _matrix[i][j]);
		}
		printf("\n");
	}
}

double& Matrix::operator()(const unsigned int &row_number, const unsigned int &column_number) {
	return _matrix[row_number][column_number];
}

unsigned int Matrix::row_size(void) {
	return _row_size;
}

unsigned int Matrix::column_size(void) {
	return _column_size;
}


Matrix Matrix::operator+(double parameter) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] + parameter;
		}
	}
	return result;
}

Matrix Matrix::operator-(double parameter) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] - parameter;
		}
	}
	return result;
}

Matrix Matrix::operator*(double parameter) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] * parameter;
		}
	}
	return result;
}

Matrix Matrix::operator/(double parameter) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] / parameter;
		}
	}
	return result;
}

Matrix Matrix::trans(void) {
	Matrix result(_column_size, _row_size);
	for (int i = 0; i < _column_size; i++) {
		for (int j = 0; j < _row_size; j++) {
			result(i, j) = _matrix[j][i];
		}
	}
	return result;
}

Matrix Matrix::operator+(Matrix &X) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] + X(i, j);
		}
	}
	return result;
}

Matrix Matrix::operator-(Matrix &X) {
	Matrix result(_row_size, _column_size);
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < _column_size; j++) {
			result(i, j) = _matrix[i][j] - X(i, j);
		}
	}
	return result;
}

Matrix Matrix::operator*(Matrix &X) {
	Matrix result(_row_size, X.column_size());
	double dummy = 0.0;
	for (int i = 0; i < _row_size; i++) {
		for (int j = 0; j < X.column_size(); j++) {
			for (int k = 0; k < _column_size; k++) {
				dummy += _matrix[i][k] * X(k, j);
			}
			result(i, j) = dummy;
			dummy = 0.0;
		}
	}
	return result;
}

Matrix Matrix::find_minor(unsigned int n, unsigned int k) {
	Matrix result(_row_size - 1, _column_size - 1);
	for (int i = 0; i < _row_size - 1; i++) {
		if (i < n) {
			for (int j = 0; j < _column_size - 1; j++) {
				if (j < k) {
					result(i, j) += _matrix[i][j];
				}
				else {
					result(i, j) += _matrix[i][j + 1];
				}
			}
		}
		else {
			for (int j = 0; j < _column_size - 1; j++) {
				if (j < k) {
					result(i, j) += _matrix[i + 1][j];
				}
				else {
					result(i, j) += _matrix[i + 1][j + 1];
				}
			}
		}
	}
	return result;
}

double Matrix::det(void) {
	if (_row_size == 1) {
		return _matrix[0][0];
	}
	else if (_row_size == 2) {
		return _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
	}
	else {
		double result = 0;
		for (int i = 0; i < _row_size; i++) {
			result += pow(-1.0, i)*_matrix[i][0] * find_minor(i, 0).det();
		}
		return result;
	}
}

Matrix Matrix::inv(void) {
	Matrix result(_row_size, _column_size);
	double determinant = det();
	if (_row_size == 1) {
		result(0, 0) = 1 / determinant;
	}
	else if (_row_size == 2) {
		result(0, 0) = _matrix[1][1] / determinant;
		result(0, 1) = _matrix[0][1] / -determinant;
		result(1, 0) = _matrix[1][0] / -determinant;
		result(1, 1) = _matrix[0][0] / determinant;
	}
	else {
		for (int i = 0; i < _row_size; i++) {
			for (int j = 0; j < _column_size; j++) {
				result(i, j) = pow(-1.0, (i + j))*find_minor(i, j).det();
			}
		}
		result = result.trans();
		result = result / determinant;
	}
	return result;
}

