#include "Matrix.h"
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <sstream>

linalg::Matrix::Matrix(int rows) :m_rows(rows) {
    
    if (m_rows <= 0) {
        throw std::runtime_error("������� ������� �����������");
    }
    m_columns = 1;
    m_ptr = new double[m_rows * m_columns];
    for (size_t i = 0; i < m_rows * m_columns; ++i) {
        m_ptr[i] = 0.0;
    }
}

linalg::Matrix::Matrix(int rows, int cols) : m_rows(rows), m_columns(cols) {
    if (m_rows <= 0 || m_columns <= 0) {
        throw std::runtime_error("������� ������� �����������");
    }
    m_ptr = new double[m_rows * m_columns];
    for (size_t i = 0; i < m_rows * m_columns; ++i) {
        m_ptr[i] = 0.0;
    }
}

linalg::Matrix::Matrix(const Matrix& object) noexcept :m_rows(object.m_rows),m_columns(object.m_columns){//����������� �����������
     m_ptr = new double[m_rows * m_columns];
     for (size_t i = 0; i < m_rows; ++i) {
         for (size_t j = 0; j < m_columns; ++j) {
             m_ptr[i * m_columns + j] = object.m_ptr[i * m_columns + j];
         }
     }
}

linalg::Matrix::Matrix(Matrix&& object) noexcept {// ����������� �����������
    std::swap(m_ptr, object.m_ptr);
    std::swap(m_rows, object.m_rows);
    std::swap(m_columns, object.m_columns);
}

linalg::Matrix::Matrix(std::initializer_list<std::initializer_list<double>> lst) {
    m_rows = lst.size();
    m_columns = lst.begin()->size();
    m_ptr = new double[m_rows * m_columns];
    size_t index = 0;
    for (const std::initializer_list<double> row : lst) {
        if (row.size() != m_columns) {
            throw std::runtime_error("������� �������� �������");
        }

        for (double element : row) {
            m_ptr[index] = element;
            index++;
        }
    }
}

void linalg::Matrix::reshape(int rows, int cols) {
    if (rows * cols != m_rows * m_columns) {
        throw std::runtime_error("������ �������� ����������� �������� ���������� ���������");
    }
    m_rows = rows;
    m_columns = cols;
}

linalg::Matrix& linalg::Matrix::operator=(const Matrix& object) {//�������� ������������ � ������������
    if (this == &object) { // �������� �� ����������������
        return *this;
    }
    if (m_rows != object.m_rows || m_columns != object.m_columns) {
        double* ptr_temp = new double[m_rows * m_columns];
        delete[] m_ptr;
        m_ptr = ptr_temp;
        m_rows = object.m_rows;
        m_columns = object.m_columns;
    }
    for (size_t i = 0; i < m_rows * m_columns; ++i) {
        m_ptr[i] = object.m_ptr[i];
    }
    return *this;
}

 linalg::Matrix& linalg::Matrix::operator=(Matrix&& object) noexcept {//�������� ������������ � ������������
     std::swap(m_ptr, object.m_ptr);
     std::swap(m_rows, object.m_rows);
     std::swap(m_columns, object.m_columns);
     return (*this);
 }
 double& linalg::Matrix::operator()(size_t row, size_t col) {// �������� ������ ������� ��� ��������� ��������
     if (row >= m_rows || col >= m_columns) {
         throw std::runtime_error("�������� ����� ��������");
     }

     return m_ptr[row * m_columns + col];
 }

double linalg::Matrix::operator()(size_t row, size_t col) const {
    if (row >= m_rows || col >= m_columns) {
        throw std::runtime_error("�������� ����� ���������");
    }

    return m_ptr[row * m_columns + col];
 }

static size_t number_of_digits(double element, std::ios_base::fmtflags flags) {//����� ��������� ������ (flags), ������� ����� ��������� ��� �������������� ����� � ������.
    std::ostringstream stream;//��������� ������ ������ std::ostringstream, ������� ������������ ����� �����, ������������ ��� ���������������� ������ ������ � ������.
    stream.flags(flags); //����� ��������������� ��������� ����� ��� ������� stream.
    //��������� ����� ������ �� ������ �������������� ����� � ������, ������� ���������� ������ ����� �������, ������ ���� � ������ ��������� ��������������.
    stream << element;//�������� element ��������� � ����� stream. ����������� ��������� �����, ������������� �����, � �������������� ����� � ������.
    return stream.str().size();//��������� �������������� ����� � ������ ����������� � ������� ������ stream.str()
    
}

std::ostream& linalg:: operator << (std::ostream& out, const Matrix& obj) {
    if (obj.empty()) // �������� ���������� ������ ������ ��������
        return out << "|empty|\n";
    size_t* max_digits = new size_t[obj.columns()]();
    size_t num_digits = 0;
    size_t digits_max = number_of_digits(obj(0, 0), out.flags());
    for (size_t j = 0; j < obj.columns(); ++j) {
        digits_max = 0;
        for (size_t i = 0; i < obj.rows(); ++i) {
            digits_max = std::max(digits_max, number_of_digits(obj(i, j), out.flags()));;
        }
        max_digits[j] = digits_max;
    }
  
    // ����� ������� � ������������� �� ��������
    for (size_t i = 0; i < obj.rows(); ++i) {
        out << '|';
        for (size_t j = 0; j < obj.columns(); ++j) {
            out << std::setw(max_digits[j]) << obj(i, j);

            if (j != obj.columns() - 1) {
                out << ' ';
            }
        }
        out << "|\n";
    }
    delete[] max_digits;
    return out;
  
}
const linalg::Matrix linalg:: operator+(const Matrix& object1, const Matrix& object2){

    if (object1.rows() != object2.rows() || object1.columns() != object2.columns() || (object1.empty() && object2.empty())) {
        throw std::runtime_error("������� ������ �������");
    }

    if (object1.empty() && !object2.empty()) {
        return object2;
    }
    if (object2.empty() && !object1.empty()) {
        return object1;
    }
    Matrix result(object1.rows(), object1.columns());

    for (size_t i = 0; i <object1.rows(); ++i) {
        for (size_t j = 0; j <object1.columns(); ++j) {
            result(i, j) = (object1)(i, j) + object2(i, j);
        }
    }

    return result;
}

linalg::Matrix& linalg::Matrix:: operator+=(const Matrix& object) {
     if (m_rows != object.m_rows || m_columns != object.m_columns) {
         throw std::runtime_error("������� ������ �������");
     }

     for (size_t i = 0; i < m_rows; ++i) {
         for (size_t j = 0; j < m_columns; ++j) {
             (*this)(i, j) += object(i, j);
         }
     }

     return *this;
 }

const linalg::Matrix linalg:: operator-(const Matrix& object1, const Matrix& object2){
    if ( (object1.rows() != object2.rows()) || (object1.columns() != object2.columns()) || (object1.empty() && object2.empty()) ) {
        throw std::runtime_error("������� ������ �������");
    }
    if (object2.empty() && !object1.empty()) {
        return object1;
    }

    Matrix result(object1.rows(), object1.columns());

    for (size_t i = 0; i < object1.rows(); ++i) {
        for (size_t j = 0; j < object1.columns(); ++j) {
            result(i, j) = (object1)(i, j) - object2(i, j);
        }
    }

    return result;
}

linalg::Matrix& linalg::Matrix:: operator-=(const Matrix& object) {
     if (m_rows != object.m_rows || m_columns != object.m_columns|| object.empty()|| (*this).empty()) {
         throw std::runtime_error("������� ������ ��������");
     }

     for (size_t i = 0; i < m_rows; ++i) {
         for (size_t j = 0; j < m_columns; ++j) {
             (*this)(i, j) -= object(i, j);
         }
     }

     return *this; 
}

const linalg::Matrix linalg:: operator*(Matrix& object1, Matrix& object2) {
    if ((object1.columns() != object2.rows()) || object1.empty()||object2.empty()) {
        throw std::runtime_error("���������� �������� �������");
    }
    Matrix result(object1.rows(), object2.columns());
    for (size_t i = 0; i < object1.rows(); ++i) {
        for (size_t j = 0; j < object2.columns(); ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < object2.columns(); ++k) {
                sum += (object1)(i, k) * object2(k, j);
            }
            result(i, j) = sum;
        }
    }

    return result;
}
const linalg::Matrix linalg::operator*(double value, const Matrix object) {
    return Matrix{ object } *= value;
}
const linalg::Matrix linalg:: operator*(const Matrix object, double value) {
    return Matrix{ object } *= value;
}

linalg::Matrix& linalg::Matrix:: operator*=(const Matrix& object) {
    if (m_columns != object.m_rows || object.empty() || (*this).empty()) {
        throw std::runtime_error("������� ������ �����������");
    }

    Matrix result(m_rows, object.m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < object.m_columns; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m_columns; ++k) {
                sum += (*this)(i, k) * object(k, j);
            }
            result(i, j) = sum;
        }
    }
    *this = result;
    return *this;
}

linalg::Matrix& linalg::Matrix:: operator*=(double value) {
    if ((*this).empty() == true) {
        throw std::runtime_error("������� ������,��� ������ �������� �� �����");
    }
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            (*this)(i, j) *= value;
        }
    }
    return *this;
}

bool linalg::operator==(const Matrix& object1, const Matrix& object2) noexcept{


   int n1 = object1.rows();
   int m1 = object1.columns();
   int n2 = object2.rows();
   int m2 = object2.columns();


    if (n1 != n2 || m1 != m2) {
        return false; // �.�. ����������� �� ���������
    }

    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < m1; ++j) {
            if ((object1)(i, j) != object2(i, j)) {
                return false; //�� ��������� ���� �� 1 �������
            }
        }
    }

    return true; 
}

bool linalg:: operator!=(const Matrix& object1, const Matrix& object2)noexcept{
    return !(object1==object2);// ������� �� (�������� ����������)=>������� ������������, ����� �� ������ ������ ���
}

double linalg::Matrix:: norm() const {
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("������� ������� �����������");
    }

    double sum = 0.0;
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            sum += m_ptr[i * m_columns + j] * m_ptr[i * m_columns + j];
        }
    }
    return sqrt(sum);
}
double linalg::Matrix:: trace() const {
    if (m_rows != m_columns){
        throw std::runtime_error("������� �� ����������, ��������� ���� ������");
    }
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("������� ������� �����������");
    }
    double tracesum = 0.0;

    for (size_t i = 0; i < m_rows; ++i) {
        tracesum += m_ptr[i * m_columns + i];
    }
    return tracesum;
}

double linalg::Matrix::determinant() const {
    if (m_rows != m_columns) {
        throw std::runtime_error("������� �� ����������");
    }

    if (m_rows == 1) {
        return m_ptr[0]; // ��� ������� 1x1 ������������ ����� � ������������� ��������
    }

    if (m_rows == 2) {
        return m_ptr[0] * m_ptr[3] - m_ptr[1] * m_ptr[2];
    }

    Matrix copy(*this);
    copy.gauss_forward();
    double det = 1.0;
    for (size_t i = 0; i < copy.m_rows; ++i)
    {
        for (size_t j = 0; j < copy.m_columns; ++j) {
            if (i == j) {
                det *= copy(i, j);
            }
        }
    }
    return det;
}

void linalg::Matrix::swapRows(size_t row_index, size_t max_row, size_t col_index) {
    // ������ ������� ������� ������ � ������ � ������������ ���������
    for (size_t i = col_index; i < m_columns; ++i) {
        double temp = (*this)(row_index, i);
        (*this)(row_index, i) = (*this)(max_row, i);
        (*this)(max_row, i) = -temp; // ������ ���� ������
    }
}
/*void linalg::Matrix::swapRows(size_t row1, size_t row2) {
    // ������ ������� ������� ������ � ������ � ������������ ���������
    for (size_t i = row1; i < m_columns; ++i) {
        double temp = (*this)(row1, i);
        (*this)(row1, i) = (*this)(row2, i);
        (*this)(row2, i) = temp; // ������ ���� ������
    }
}*/

/*void linalg::Matrix::gauss_forward() {
    if ((*this).empty()) {
        throw std::runtime_error("������� �����, �� ��������� ��������� ������ ��� ������.");
    }
    for (size_t nowrow = 0; nowrow < m_rows; ++nowrow) {
        // ����� ������ ��������� ������� � ������� �������
        size_t nozerorow = nowrow;//nozerow �������� �� ������, ����������� �� ��������

        while (nozerorow < m_rows && (*this)(nozerorow, nowrow) == 0.0) {
            ++nozerorow;
        }

        if (nozerorow == m_rows) {// �� ������� ��������� ��������� � �������=> ��������� � ���������� �������
            continue;
        }

        if (nozerorow != nowrow) {
            swapRows(nozerorow, nowrow);// �������� ������� ������, ����� ��������� ������ � ��������� ���������  �� ���������
        }
        for (size_t j = nowrow + 1; j < m_rows; ++j) {
            double divide = (*this)(j, nowrow) / (*this)(nowrow, nowrow);
            for (size_t i = nowrow; i < m_columns; ++i) {//�������� �� ��������(��� �� ��������� ����� ������)
                (*this)(j, i) -= divide * (*this)(nowrow, i);
            }
        }
    }
}*/

 void linalg::Matrix::gauss_forward() {
        const double mistake = 1e-6;
        if ((*this).empty()) {
            throw std::runtime_error("������� ������. ���������� ��������� ������ ��� ������ ������");
        }
        size_t row_index = 0, col_index = 0;

        while (row_index < m_rows && col_index < m_columns) {
            // ������� ������������ ������� � ������� �������
            size_t max_row = row_index;
            double max_element = fabs((*this)(row_index, col_index));

            for (size_t i = row_index; i < m_rows; ++i) {
                double current_element = fabs((*this)(i, col_index));
                if (current_element > max_element) {
                    max_element = current_element;
                    max_row = i;
                }
            }

            if (max_element <=mistake) {
                for (size_t i = row_index; i < m_rows; ++i) {
                    (*this)(i, col_index) = 0.0;
                }
                ++col_index; // ��������� � ���������� �������
            }
            else {
                if (max_row != row_index) {
                    // ������ ������� ������� ������ � ������ � ������������ ���������
                    swapRows(row_index, max_row, col_index);
                }

                for (size_t i = row_index + 1; i < m_rows; ++i) {
                    double factor = -((*this)(i, col_index) / (*this)(row_index, col_index));
                    (*this)(i, col_index) = 0.0;

                    for (size_t j = col_index + 1; j < m_columns; ++j) {
                        (*this)(i, j) += factor * (*this)(row_index, j);
                    }
                }

                ++row_index;
                ++col_index;
            }
        }

}

void linalg::Matrix::gauss_backward() {
    if (empty()) {
        throw std::runtime_error("������� ������. ���������� ��������� �������� ��� ������ ������");
    }

    for (int row = m_rows - 1; row >= 0; --row) {
        int nowColumn = -1; // �������������� ������� ��� -1

        // ������� ������ ��������� ������� � ������� ������
        for (size_t col = 0; col < m_columns; ++col) {
            if ((*this)(row, col)!=0) {
                nowColumn = col;
                break;
            }
        }

        if (nowColumn == -1) {
            continue;
        }

        // ����������� ������, ����� ������� ������� ���� ����� 1
        double pivotValue = (*this)(row, nowColumn);
        for (size_t col = 0; col < m_columns; ++col) {
                (*this)(row, col) /= pivotValue;
        }

        // �������� �������� ���� �������� �������� � �������
        for (int upper = row - 1; upper >= 0; --upper) {
            double value = - (*this)(upper, nowColumn);
            for (size_t col = 0; col < m_columns; ++col) {
                if (col == nowColumn) {
                    (*this)(upper, col) = 0;
                }
                else {
                    (*this)(upper, col) += value * (*this)(row, col);
                }
            }
        }
    }
}

int linalg::Matrix::rank() const {
    if (m_rows == 0 || m_columns == 0 || (*this).empty()) {
        throw std::invalid_argument("������� ������� �����������");
    }
    int rank = 0;
    Matrix temp(*this);//�����,������ ��� ����� ������ ������ �������
    temp.gauss_forward();
    for (size_t i = 0; i < temp.m_rows; ++i) {
        bool nonullrow = false;
        for (size_t j = 0; j < temp.m_columns; ++j) {
            if (temp(i, j) != 0.0) {
                nonullrow = true;
                break;
            }
        }
        if (nonullrow) {
            ++rank;
        }
    }
    return rank;
}

linalg:: Matrix linalg::concatenate(const Matrix &left, const Matrix &right) {

    if (left.empty()||right.empty()) {
        throw std::runtime_error("������� ������. ���������� ��������� �������� ��� ������ ������");
    }
    if (left.rows() != right.rows()) {
        throw std::runtime_error("���������� ���������� ������� � ������ ����������� �����.");
    }

    int commonCols = left.columns() + right.columns();
    Matrix result(left.rows(), commonCols);

    for (int i = 0; i < left.rows(); ++i) {
        for (int j = 0; j < left.columns(); ++j) {
            result(i, j) = left(i, j);
        }
        for (int j = 0; j < right.columns(); ++j) {
            result(i, left.columns() + j) = right(i, j);
        }
    }

    return result;
}
linalg::Matrix linalg::transpose(Matrix matr) {
    if (matr.empty()) {
        throw std::runtime_error("���������� ���������������");
    }
    int n = matr.rows();
    int m= matr.columns();
    Matrix result(m, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            result(j, i) = matr(i, j);
        }
    }

    return result;
}

linalg::Matrix linalg:: invert(Matrix matr) {
    const double mistake = 1e-6;
    if (matr.empty()) {
        throw std::runtime_error("������� ������");
    }
    if (matr.columns() != matr.rows()){
        throw std::runtime_error("������� �� ����������, ���������� ����� ��������");
    }
    if (matr.determinant() == 0) {
        throw std::runtime_error("���������� ����� �������� �������");
    }

    size_t n = matr.rows();
    size_t m = matr.columns();
    Matrix backMatr(n, 2 * n); 


    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < 2 * n; ++j) {
            if (j < n) {
                backMatr(i, j) = matr(i, j);
            }
            else if (j == (i + n)) {
                backMatr(i, j) = 1.0;
            }
            else {
                backMatr(i, j) = 0.0;
            }
        }
    }

    backMatr.gauss_forward();
    backMatr.gauss_backward();




    Matrix result(n, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result(i, j) = backMatr(i, j + n);
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(result(i, j)) <= mistake) {
                result(i, j) = 0;
            }
        }
    }

    return result;
}

linalg::Matrix linalg::power(Matrix& matr, int pow) {
    if (matr.empty()) {
        throw std::runtime_error("������� ������");
    }
    if (matr.columns() != matr.rows()) {
        throw std::runtime_error("������� �� ����������, ���������� �������� � �������");
    }
    Matrix MatrPow(matr.columns(), matr.rows());
    int n = matr.columns();
    if (pow == 0) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    MatrPow(i, j) = 1.0;
                }
                else {
                    MatrPow(i, j) = 0.0;
                }
    
            }
        }
        return MatrPow;
    }

    if (pow == 1) {
        return matr;
    }

    if (pow > 0) {
        Matrix MatrPow = matr;
        for (int i = 0; i < pow - 1; ++i) {
            MatrPow = MatrPow * matr;
        }
        return MatrPow;
    }

    if (pow < 0) {
        MatrPow = invert(matr);
        for (int i = 0; i < -pow - 1; ++i) {
            MatrPow = MatrPow * matr;
        }
        return MatrPow;

    }
}

linalg::Matrix linalg::solve(const Matrix& matr_a, const Matrix& vec_f) {
    
    Matrix solution = concatenate(matr_a, vec_f);
    std::cout << solution.rank();
    std::cout << matr_a.rank();
    if (matr_a.empty() || vec_f.empty()){
        throw std::runtime_error("������ ����� ������");
    }
    if (solution.rank() != matr_a.rank()){
        throw std::runtime_error("������� �� ����� �������");
    }
    if ((solution.rank() == matr_a.rank()) && (solution.rank()!=matr_a.columns())){
        throw std::runtime_error("������� ����� ���������� ����� �������");
    }
    solution.gauss_forward();
    solution.gauss_backward();
    linalg::Matrix solution_vector(vec_f.rows(), 1);
    for (size_t i = 0; i < solution_vector.rows(); i++) {
        solution_vector(i, 0) = solution(i, solution.columns() - 1);
    }

    return solution_vector;
}

/*�������������*/

linalg::Matrix linalg::Matrix::operator-() const { //�������� �������� ������
    linalg::Matrix result(m_rows, m_columns);
    for (size_t i = 0; i < m_rows; ++i) {
        for (size_t j = 0; j < m_columns; ++j) {
            result(i, j) = -(*this)(i, j);
        }
    }
    return result;
}

double linalg::Matrix::getMinor(int firstRow, int firstCol, int n) const {
    if (firstRow < 0 || firstCol < 0 || firstRow + n > m_rows || firstCol + n > m_columns) {
        throw std::runtime_error("��������� �����������");
    }

    linalg::Matrix minorMatrix(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            minorMatrix(i, j) = (*this)(firstRow + i, firstCol + j);
        }
    }
    return minorMatrix.determinant();
}
