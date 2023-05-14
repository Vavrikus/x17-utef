// C++ dependencies
#include <vector>
#include <iostream>

// X17 dependencies
#include "Matrix.h"

namespace X17
{   
    //// Public methods.

    template <int M, int N>    
    inline const double& Matrix<M,N>::at(int row, int column) const
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

    template <int M, int N>
    inline double& Matrix<M,N>::at(int row, int column)
    {
        #ifdef DEBUG
        if (row < 0 || row >= M || column < 0 || column >= N)
        {
            std::cerr << "ERROR: Invalid matrix element (" << row << "," << column << ") of " << M << "x" << N << " matrix.\n";
            return elements[0];
        }
        #endif
        return elements[row * N + column];
    }

    template <int M, int N>
    std::vector<double> Matrix<M,N>::GetColumn(int c) const
    {
        // Check if the column index is valid.
        if(c >= 0 && c < N)
        {
            std::vector<double> column(M);

            for (int r = 0; r < M; r++) column[r] = this->at(r, c);
            return column;
        }
        else
        {
            std::cerr << "ERROR: Invalid matrix column index.\n";
            return std::vector<double>();
        }
    }

    template <int M, int N>
    void Matrix<M,N>::Print() const
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

    template <int M, int N>
    void Matrix<M,N>::SwapRows(int row1, int row2)
    {
        // Check if the row indexes are valid.
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

    template <int M, int N>
    void Matrix<M,N>::Reduce() 
    {
        for (int c = 0; c < M; c++) 
        {
            // Find row with non-zero c-th element.
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

            // Normalize selected row.
            double first = this->at(c, c);
            if (first == 0) continue; // Skip division by zero.
            for (int c2 = c; c2 < N; c2++) this->at(c, c2) /= first;

            // Subtracting rows.
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

    template <int M, int N>
    Matrix<M,N> operator+(const Matrix<M,N>& A, const Matrix<M,N>& B)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M * N; i++) 
        {
            result.elements[i] = A.elements[i] + B.elements[i];
        }

        return result;
    }

    template <int M, int N>
    Matrix<M,N> operator*(double d, const Matrix<M,N>& A)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M * N; i++)
        {
            result.elements[i] = d * A.elements[i];
        }
        return result;
    }

    template <int M, int N>
    Matrix<M,N> operator/(const Matrix<M,N>& A, double d)
    {
        Matrix<M,N> result;
        for (int i = 0; i < M * N; i++)
        {
            result.elements[i] = A.elements[i] / d;
        }
        return result;
    }

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