#pragma once

// C++ dependencies
#include <vector>

namespace X17
{
    /// @brief A templated Matrix class with M rows and N columns. This class represents a matrix of size M x N,
    ///        where the matrix elements are stored in a contiguous array of size M*N. The class provides methods
    ///        for performing basic mathematical operations on matrices, such as addition and multiplication by a scalar.
    /// @tparam M The number of rows in the matrix.
    /// @tparam N The number of columns in the matrix.     
    template <int M, int N>
    struct Matrix
    {
        double elements[M*N]; // The array storing the matrix elements in row-major order.

        /// @brief Default constructor.
        Matrix() = default;

        /// @brief Constructor that initializes the matrix from an array of size M*N.
        /// @param arr The array of size M*N containing the matrix elements in row-major order.
        Matrix(const double (&arr)[M*N])
        {
            copy(std::begin(arr), std::end(arr), std::begin(elements));
        }

        /// @brief Addition operator that adds a matrix to the current matrix.
        /// @param A The matrix to add.
        void operator+=(const Matrix<M,N>& A)
        {
            for (int i = 0; i < M*N; i++) this->elements[i] += A.elements[i];
        }

        /// @brief Multiply each element in the matrix by a scalar.
        /// @param d The scalar to multiply the elements by.
        void operator*=(const double& d)
        {
            for (double& element : elements) element *= d;
        }
    
        /// Access the matrix element at the specified row and column.
        /// @param row The row of the element to access, starting from 0.
        /// @param column The column of the element to access, starting from 0.
        /// @return A constant reference to the matrix element at the specified row and column.
        /// If the specified row or column is out of bounds, returns a reference to the first element and prints an error message.
        /// @note If DEBUG is defined, this method checks if the specified row and column are within bounds and prints an error message if not.
        /// If DEBUG is not defined, no bounds checking is performed for efficiency.
        const double& at(int row, int column) const;

        /// Access the matrix element at the specified row and column.
        /// @param row The row of the element to access, starting from 0.
        /// @param column The column of the element to access, starting from 0.
        /// @return A reference to the matrix element at the specified row and column.
        /// If the specified row or column is out of bounds, returns a reference to the first element and prints an error message.
        /// @note If DEBUG is defined, this method checks if the specified row and column are within bounds and prints an error message if not.
        /// If DEBUG is not defined, no bounds checking is performed for efficiency.
        double& at(int row, int column);

        /// @brief Returns a vector representing the specified column of the matrix.
        /// @param c The index of the column to retrieve.
        /// @return A vector representing the specified column of the matrix, or an empty vector if the index is out of bounds.
        std::vector<double> GetColumn(int c) const;

        /// @brief Prints the matrix to the console.
        void Print() const;

        /// @brief Swaps the elements of two rows in the matrix.
        /// @param row1 The index of the first row to swap.
        /// @param row2 The index of the second row to swap.
        void SwapRows(int row1, int row2);

        /// @brief Reduces the matrix to row echelon form using Gaussian elimination.
        void Reduce();
    };

    /// @brief Adds two matrices element-wise and returns the result.
    /// @tparam M the number of rows in the matrices.
    /// @tparam N the number of columns in the matrices.
    /// @param A The first matrix.
    /// @param B The second matrix.
    /// @return The sum of the two matrices.
    template <int M, int N>
    Matrix<M,N> operator+(const Matrix<M,N>& A, const Matrix<M,N>& B);

    /// @brief Multiplies a scalar with a matrix element-wise.
    /// @tparam M The number of rows in the matrix.
    /// @tparam N The number of columns in the matrix.
    /// @param d The scalar to multiply with.
    /// @param A The matrix to multiply with the scalar.
    /// @return The resulting matrix of the element-wise multiplication.
    template <int M, int N>
    Matrix<M,N> operator*(const double& d, const Matrix<M,N>& A);

    /// @brief Divide each element of a matrix by a scalar.
    /// @tparam M Number of rows of the matrix.
    /// @tparam N Number of columns of the matrix.
    /// @param A Matrix to be divided.
    /// @param d Scalar to divide the matrix by.
    /// @return Matrix result of dividing each element of A by d.
    template <int M, int N>
    Matrix<M,N> operator/(const Matrix<M,N>& A,const double& d);

    /// @brief Multiplies two matrices A and B, and returns the result C = A*B.
    /// @tparam M Number of rows of matrix A and C.
    /// @tparam N Number of columns of matrix A and rows of matrix B.
    /// @tparam P Number of columns of matrix B and C.
    /// @param A First matrix to multiply.
    /// @param B Second matrix to multiply.
    /// @return Matrix<M,P> Resulting matrix of size MxP.
    template<int M,int N,int P>
    Matrix<M,P> operator*(const Matrix<M,N>& A, const Matrix<N,P>& B);
} // namespace X17