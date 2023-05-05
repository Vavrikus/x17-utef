#pragma once

#include <iostream>
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
        const double& at(int row, int column) const
        {
            #ifdef DEBUG
            if (row < 0 || row >= M || column < 0 || column >= N)
            {
                std::cerr << "ERROR: Invalid matrix element (" << row << "," << column << ") of " << M << "x" << N << " matrix.\n";
                return elements[0];
            }
            #endif
            return elements[row*N+column];
        }

        /// Access the matrix element at the specified row and column.
        /// @param row The row of the element to access, starting from 0.
        /// @param column The column of the element to access, starting from 0.
        /// @return A reference to the matrix element at the specified row and column.
        /// If the specified row or column is out of bounds, returns a reference to the first element and prints an error message.
        /// @note If DEBUG is defined, this method checks if the specified row and column are within bounds and prints an error message if not.
        /// If DEBUG is not defined, no bounds checking is performed for efficiency.
        double& at(int row, int column)
        {
            #ifdef DEBUG
            if (row < 0 || row >= M || column < 0 || column >= N)
            {
                std::cerr << "ERROR: Invalid matrix element (" << row << "," << column << ") of " << M << "x" << N << " matrix.\n";
                return elements[0];
            }
            #endif
            return elements[row*N+column];
        }

        /// @brief Returns a vector representing the specified column of the matrix.
        /// @param c The index of the column to retrieve.
        /// @return A vector representing the specified column of the matrix, or an empty vector if the index is out of bounds.
        std::vector<double> GetColumn(int c) const
        {
            // Check if the column index is valid
            if(c >= 0 && c < N)
            {
                std::vector<double> column(M);

                for (int r = 0; r < M; r++) column[r] = this->at(r, c);
                return column;
            }
            else
            {
                std::cerr << "ERROR: Invalid matrix column index.\n";
                return vector<double>();
            }
        }

        /// @brief Prints the matrix to the console.
        void Print() const
        {
            for (int r = 0; r < M; r++)
            {
                for (int c = 0; c < N; c++)
                {
                    std::cout << this->at(r,c) << " ";
                }
                std::cout << "\n";
            }        
        }

        /// @brief Swaps the elements of two rows in the matrix.
        /// @param row1 The index of the first row to swap.
        /// @param row2 The index of the second row to swap.
        void SwapRows(int row1, int row2)
        {
            // Check if the row indexes are valid
            if (row1 > -1 && row1 < M && row2 > -1 && row2 < M)
            {
                for (int c = 0; c < N; c++)
                {
                    double temp = this->at(row1,c);
                    this->at(row1,c) = this->at(row2,c);
                    this->at(row2,c) = temp;
                }            
            }

            else std::cerr << "ERROR: Invalid matrix row number.\n";     
        }

        /// @brief Reduces the matrix to row echelon form using Gaussian elimination.
        void Reduce() 
        {
            for (int c = 0; c < M; c++) 
            {
                // Find row with non-zero c-th element
                bool not_zero = false;
                for (int r = c; r < M; r++) 
                {
                    if (this->at(r, c) != 0) 
                    {
                        not_zero = true;
                        if (r != c) this->SwapRows(r, c);
                        break;
                    }
                }
                if (!not_zero) std::cerr << "WARNING: Singular matrix.\n";

                // Normalize selected row
                double first = this->at(c, c);
                if (first == 0) continue; // Skip division by zero
                for (int c2 = c; c2 < N; c2++) this->at(c, c2) /= first;

                // Subtracting rows
                for (int r = 0; r < M; r++) 
                {
                    if (r != c) 
                    {
                        double factor = this->at(r, c);
                        for (int c2 = c; c2 < N; c2++) {
                            this->at(r, c2) -= factor * this->at(c, c2);
                        }
                    }
                }
            }
        }
    };

    /// @brief Adds two matrices element-wise and returns the result.
    /// @tparam M the number of rows in the matrices.
    /// @tparam N the number of columns in the matrices.
    /// @param A The first matrix.
    /// @param B The second matrix.
    /// @return The sum of the two matrices.
    template <int M, int N>
    Matrix<M,N> operator+(const Matrix<M,N>& A, const Matrix<M,N>& B)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M*N; i++) 
        {
            result.elements[i] = A.elements[i] + B.elements[i];
        }

        return result;
    }

    /// @brief Multiplies a scalar with a matrix element-wise.
    /// @tparam M The number of rows in the matrix.
    /// @tparam N The number of columns in the matrix.
    /// @param d The scalar to multiply with.
    /// @param A The matrix to multiply with the scalar.
    /// @return The resulting matrix of the element-wise multiplication.
    template <int M, int N>
    Matrix<M,N> operator*(const double& d, const Matrix<M,N>& A)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M*N; i++)
        {
            result.elements[i] = d * A.elements[i];
        }
        return result;
    }

    /// @brief Divide each element of a matrix by a scalar.
    /// @tparam M Number of rows of the matrix.
    /// @tparam N Number of columns of the matrix.
    /// @param A Matrix to be divided.
    /// @param d Scalar to divide the matrix by.
    /// @return Matrix result of dividing each element of A by d.
    template <int M, int N>
    Matrix<M,N> operator/(const Matrix<M,N>& A,const double& d)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M*N; i++)
        {
            result.elements[i] = A.elements[i] / d;
        }
        return result;
    }

    /// @brief Multiplies two matrices A and B, and returns the result C = A*B.
    /// @tparam M Number of rows of matrix A and C.
    /// @tparam N Number of columns of matrix A and rows of matrix B.
    /// @tparam P Number of columns of matrix B and C.
    /// @param A First matrix to multiply.
    /// @param B Second matrix to multiply.
    /// @return Matrix<M,P> Resulting matrix of size MxP.
    template<int M,int N,int P>
    Matrix<M,P> operator*(const Matrix<M,N>& A, const Matrix<N,P>& B)
    {
        Matrix<M,P> result;
        for (int r = 0; r < M; r++)
        {
            for (int c = 0; c < P; c++)
            {
                double sum = 0;
                for (int i = 0; i < N; i++) sum += A.at(r,i)*B.at(i,c);
                result.at(r,c) = sum;
            }        
        }
        
        return result;
    }
} // namespace X17