#pragma once
#include <initializer_list>
#include <iostream>
namespace linalg {
	class Matrix {
	public:
		Matrix() : m_ptr(nullptr), m_rows(0), m_columns(0) {}//дефолтный конструктор
		Matrix(int rows, int cols);//конструктор заполняет нулями
		Matrix(int rows);
		Matrix(std::initializer_list<std::initializer_list<double>> lst);
		Matrix(const Matrix& object)noexcept;//конструктор копирования
		Matrix(Matrix&& object)noexcept;
		int rows() const noexcept { return m_rows; }
		int columns() const noexcept { return m_columns; }
		bool empty() const { return m_ptr == nullptr; }
		void reshape(int rows, int cols);
		double norm() const;
		double trace() const;
		void gauss_forward();
		void gauss_backward();
		double determinant() const;
		int rank() const;
		double getMinor(int firstRow, int firstCol, int n) const;
		Matrix& operator=(const Matrix& object);
		Matrix& operator=(Matrix&& object)noexcept;
		double& operator()(size_t row, size_t col);
		double operator()(size_t row, size_t col) const;
		Matrix& operator +=(const Matrix& object);
		Matrix& operator -=(const Matrix& object);
		Matrix& operator*=(const Matrix& object);
		Matrix& operator*=(double digit);
		~Matrix() { delete[] m_ptr; }
		Matrix operator -() const;
	private:
		void swapRows(size_t row_index, size_t max_row, size_t col_index);
		/*void swapRows(size_t row1, size_t row2);*/
		double* m_ptr=nullptr;
		int m_rows=0;
		int m_columns=0;

	};
	const Matrix operator*(const Matrix object, double value);
	const Matrix operator*(double value, const Matrix object);
	const Matrix operator*(Matrix& object2, Matrix& object1);
	const Matrix operator+(const Matrix& object1, const Matrix& object2);
	const Matrix operator-(const Matrix& object1, const Matrix& object2);
	bool operator==(const Matrix& object1, const Matrix& object2) noexcept;
	bool operator!=(const Matrix& object1, const Matrix& object2) noexcept;
	std::ostream& operator<<(std::ostream& out, const Matrix& matrix);
	Matrix concatenate(const Matrix& left, const Matrix& right);
	Matrix transpose(Matrix matr);
	Matrix invert(Matrix matr);
	Matrix power(Matrix& matr, int pow);
	Matrix solve(const Matrix& matr_a, const Matrix& vek_f);
}