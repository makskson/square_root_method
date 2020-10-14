#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

const int SIZE = 5;

double A[SIZE][SIZE] = {
	{0.8894, 0.0000, -0.2323, 0.1634, 0.2723},
	{-0.0545, 0.5808, 0.0000, -0.1107, 0.0363},
	{0.0182, -0.1634, 1.0527, 0.0200, 0.0635},
	{0.0545, 0.0000, -0.1325, 1.0527, 0.0000},
	{0.0363, -0.0545, 0.2632, -0.0218, 0.7623}
};
double b[SIZE] = {
	4.2326, -4.1037, -2.6935, 1.6916, 3.1908
};
double Atr[SIZE][SIZE];
double Anew[SIZE][SIZE];
double bnew[SIZE];
double y[SIZE];
double S[SIZE][SIZE];
double Str[SIZE][SIZE];
double x[SIZE];
double n[SIZE];

void transpose(double matrix[SIZE][SIZE], double newMatrix[SIZE][SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			newMatrix[j][i] = A[i][j];
		}
	}
}
void separate() {
	cout << "__________________________________________________________________\n";
}

void show_matrix(double matrix[SIZE][SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			cout << setw(13) << matrix[i][j];
		}
		cout << endl;
	}
	separate();
}
void show_advanced_matrix(double matrix[SIZE][SIZE], double vector[SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			cout << setw(13) << matrix[i][j];
		}
		cout << "|" << setw(13) << vector[i];
		cout << endl;
	}
	separate();
}
void show_vector(double vector[SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		cout << "[" << i + 1 << "] " << vector[i] << endl;
	}
	separate();
}

void find_solution(double matrix[SIZE][SIZE], double b[SIZE]) {
	// Нахождение матрицы S.
	S[0][0] = sqrt(matrix[0][0]);
	for (int i = 1; i < SIZE; i++) {
		S[0][i] = matrix[0][i] / S[0][0];
	}
	for (int i = 1; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			if (i == j) {
				double sum = 0;
				for (int k = 0; k <= i - 1; k++) {
					sum += S[k][i] * S[k][i];
				}
				S[i][j] = sqrt(matrix[i][j] - sum);
			}
			else if (i > j) {
				S[i][j] = 0;
			}
			else {
				double sum = 0;
				for (int k = 0; k <= i - 1; k++) {
					sum += S[k][i] * S[k][j];
				}
				S[i][j] = (matrix[i][j] - sum) / S[i][i];
			}
		}
	}
	transpose(S, Str);
	// Нахождение вектора у.
	y[0] = b[0] / S[0][0];
	for (int i = 1; i < SIZE; i++) {
		double sum = 0;
		for (int k = 0; k <= i - 1; k++) {
			sum += S[k][i] * y[k];
		}
		y[i] = (b[i] - sum) / S[i][i];
	}
	// Нахождение вектора х.
	x[SIZE - 1] = y[SIZE - 1] / S[SIZE - 1][SIZE - 1];
	for (int i = SIZE - 2; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k < SIZE; k++) {
			sum += S[i][k] * x[k];
		}
		x[i] = (y[i] - sum) / S[i][i];
	}
}

void multiply_matrix_vector(double matrix[SIZE][SIZE], double vector[SIZE], double answer[SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		double sum = 0;
		for (int j = 0; j < SIZE; j++) {
			sum += matrix[i][j] * vector[j];
		}
		answer[i] = sum;
	}
}
void multiply_matrix_matrix(double matrix1[SIZE][SIZE], double matrix2[SIZE][SIZE], double answer[SIZE][SIZE]) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			double sum = 0;
			for (int k = 0; k < SIZE; k++) {
				sum += matrix1[i][k] * matrix2[k][j];
			}
			answer[i][j] = sum;
		}
	}
}

int main() {
	setlocale(LC_ALL, "rus");
	cout.precision(4);
	cout << scientific;
	cout.fixed;

	cout << "Первоначальная матрица А:\n";
	show_matrix(A);
	cout << "Первоначальный вектор b:\n";
	show_vector(b);
	transpose(A, Atr);
	multiply_matrix_matrix(Atr, A, Anew);
	cout << "Симметричекая матрица А:\n";
	show_matrix(Anew);
	multiply_matrix_vector(Atr, b, bnew);
	cout << "Новый вектор b:\n";
	show_vector(bnew);
	find_solution(Anew, bnew);
	cout << "Вектор х:\n";
	show_vector(x);
	cout << "Вектор у:\n";
	show_vector(y);
	cout << "Матрица S:\n";
	show_matrix(S);
	double det = 1.;
	for (int i = 0; i < SIZE; i++) {
		det *= S[i][i];
	}
	multiply_matrix_vector(Anew, x, n);
	for (int i = 0; i < SIZE; i++) {
		n[i] -= bnew[i];
	}
	cout << "Невязка:\n";
	show_vector(n);
	cout << "Определитель матрицы А:\n" << det << endl;
}