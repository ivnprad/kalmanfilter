#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <initializer_list>
#include <stddef.h>
#include <algorithm>
#include <cassert>
#include "Utils.hpp"

//#include "quaternion.h"

namespace robotics{
    //class Quaternion,
    template<typename T, size_t ROWS, size_t COLS>
    class Matrix{
        public:
            Matrix()=default;
            explicit Matrix(const T& value) { // this is not recognized not sure why!
                std::fill_n(&_data[0][0], ROWS*COLS, value);
            }

            //#TODO MISSING DEEP COPY CONSTRUCTOR AND DEEP COPY ASSIGNMENT CONSTRUCTOR
            Matrix(std::initializer_list<std::initializer_list<T>> init); 

            //RUN-TIME
            T& at(size_t row, size_t col);
            const T& at(size_t row, size_t col) const;

            //COMPILE-TIME
            constexpr T& atCT(size_t row, size_t col) {
                return _data[row][col];
            }
            constexpr const T& atCT(size_t row, size_t col) const {
                return _data[row][col];
            }            

            //comparators
            bool operator==(const Matrix<T,ROWS,COLS>& other) const;
            bool operator!=(const Matrix<T,ROWS,COLS>& other) const;


            //compound assigment operators
            Matrix<T,ROWS,COLS>& operator+=(const Matrix<T,ROWS,COLS>& other);
            Matrix<T,ROWS,COLS>& operator-=(const Matrix<T,ROWS,COLS>& other);
            Matrix<T,ROWS,COLS>& operator/=(const Matrix<T,ROWS,COLS>& other);
            Matrix<T,ROWS,COLS>& operator*=(T scalar);
            Matrix<T,ROWS,COLS>& operator/=(T scalar);
            Matrix<T,ROWS,COLS>& operator+=(T scalar);
            Matrix<T,ROWS,COLS>& operator-=(T scalar);

            //binary arithmetic operators
            Matrix<T,ROWS,COLS> operator*(T scalar) const;
            Matrix<T,ROWS,COLS> operator+(T scalar) const;
            Matrix<T,ROWS,COLS> operator-(T scalar) const;
            Matrix<T,ROWS,COLS> operator+(const Matrix<T,ROWS,COLS>& other)const;
            Matrix<T,ROWS,COLS> operator-(const Matrix<T,ROWS,COLS>& other)const;
            Matrix<T,ROWS,COLS> operator/(const Matrix<T,ROWS,COLS>& other)const; //#TODO this might be wrong implemented

            // all elements sum
            T sum() const;

            //
            size_t getRows() const {return ROWS;};
            //size_t getRows() {return ROWS;};
            size_t getColumns() {return COLS;};

            Matrix<T,1,COLS> getRow(const size_t row_idx);
            Matrix<T,ROWS,1> getCol(const size_t col_idx);

            // Linear Algebra
            Matrix<T,ROWS,COLS> operator*(const Matrix<T,ROWS,COLS>& other)const;
            Matrix<T,ROWS-1,COLS-1> getMinor(size_t row_idx, size_t col_idx)const;
            static Matrix<T,ROWS,COLS> identity();
            //static const Matrix<T,ROWS,COLS> GetRotMatrix(const Quaternion& q_rot); // later
            T getDeterminant();
            Matrix<T,COLS,ROWS> getTransposed() const;
            Matrix<T,COLS,ROWS> getAdjoint() const;
            Matrix<T,ROWS,COLS> getAdjugate() const;

            template <typename S = double>
            Matrix<S,ROWS,COLS> getInverse();

        private:
            T _data[ROWS][COLS]; // allocated in the stack no need for deep copy constructor this is only needed for heap pointers new delete

    };

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>::Matrix(std::initializer_list<std::initializer_list<T>> init){
        size_t row = 0;
        for(const auto& row_list: init){
            //iterate over each element of the row
            size_t col = 0;
            for (const auto& value: row_list){
                at(row,col++)=value;
            }
            ++row;
        }
    }

    // retrieve entry of Matrix at row and col 
    template<typename T, size_t ROWS, size_t COLS>
    T& Matrix<T,ROWS,COLS>::at(size_t row, size_t col){
        if (row>=ROWS){
            throw std::out_of_range("row index out of range");
        }
        if ( col >= COLS){
            throw std::out_of_range("col index out of range");
        }
        return _data[row][col];
    }

    // const retrieve entry of Matrix at row and col 
    template<typename T, size_t ROWS, size_t COLS>
    const T& Matrix<T,ROWS,COLS>::at(size_t row, size_t col)const{
        if (row >=ROWS){
            throw std::out_of_range("row index out of range");
        }
        if ( col >= COLS){
            throw std::out_of_range("col index out of range");
        }
        return _data[row][col];
    }


    //bool operator
    template<typename T, size_t ROWS, size_t COLS>
    bool Matrix<T,ROWS,COLS>::operator==(const Matrix<T,ROWS,COLS>& other) const{
        if(this==&other){// comparing with itself
            return true;
        }

        if constexpr(ROWS==0 || COLS ==0){
            return true;
        }

        return std::equal(&_data[0][0],&_data[ROWS-1][COLS-1]+1,&other._data[0][0]);
    }

    template<typename T, size_t ROWS, size_t COLS>
    bool Matrix<T,ROWS,COLS>::operator!=(const Matrix<T,ROWS,COLS>& other) const{
        return !(*this==other);
    }

    //compound assignment operators
    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator+=(const Matrix<T,ROWS,COLS>& other) {
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                this->at(i,j) += other.at(i,j);
            }
        }
        return *this;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator-=(const Matrix<T,ROWS,COLS>& other) {
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                this->at(i,j) -= other.at(i,j);
            }
        }
        return *this;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator/=(const Matrix<T,ROWS,COLS>& other) {
        for (size_t i = 0; i < ROWS; i++) {
            for (size_t j = 0; j < COLS; j++) {
            this->at(i,j) /= other.at(i,j);
            }
        }
        return *this;
    }


    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator*=(T scalar) {
        // The +1 in &(this->_data[ROWS-1][COLS-1]+1) is used to obtain a pointer one element past the end of the array. The reason this is done is 
        // because the std::for_each algorithm works on a range that includes the start element but excludes
        // the end element, so we need to add 1 to the pointer to make sure the last element is included in the range
        std::for_each(&(this->_data[0][0]),&(this->_data[ROWS-1][COLS-1])+1,[=](auto& element){
            element *= scalar;
        });
        return *this;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator+=(T scalar){
        std::for_each(&(this->_data[0][0]),&(this->_data[ROWS-1][COLS-1])+1,[=](auto& element){
            element += scalar;
        });
        return *this;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator-=(T scalar){
        std::for_each(&(this->_data[0][0]),&(this->_data[ROWS-1][COLS-1])+1,[=](auto& element){
            element -= scalar;
        });
        return *this;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS>& Matrix<T,ROWS,COLS>::operator/=(T scalar){
        assert(static_cast<int>(scalar)!=0);
        T new_scalar = 1/scalar;
        (*this)*=new_scalar;
        return *this;
    }


    // Binary arithmetic operators 
    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator*(T scalar)const{
        return Matrix<T,ROWS,COLS>(*this)*=scalar;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator+(T scalar) const{
        return Matrix<T,ROWS,COLS>(*this)+=scalar;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator-(T scalar) const{
        return Matrix<T,ROWS,COLS>(*this)-=scalar;
    }    

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator+(const Matrix<T,ROWS,COLS>& other) const{
        return Matrix<T,ROWS,COLS>(*this)+=other;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator-(const Matrix<T,ROWS,COLS>& other) const{
        return Matrix<T,ROWS,COLS>(*this)-=other;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator/(const Matrix<T,ROWS,COLS>& other)const{
        return Matrix<T,ROWS,COLS>(*this)/=other;
    }

    template<typename T, size_t ROWS, size_t COLS>
    T Matrix<T,ROWS,COLS>::sum() const{
        T ret=static_cast<T>(0);
        std::for_each(&(this->_data[0][0]),&(this->_data[ROWS-1][COLS-1])+1,[&](auto& element){
            ret+=element;
        });
        return ret;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,1,COLS> Matrix<T,ROWS,COLS>::getRow(const size_t row_idx){
        Matrix<T,1,COLS> row_vector_at_idx;
        for (size_t col_idx =0; col_idx<COLS;++col_idx){
            row_vector_at_idx.at(0,col_idx) = at(row_idx,col_idx);
        } 
        return row_vector_at_idx;   
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,1> Matrix<T,ROWS,COLS>::getCol(const size_t col_idx){
        Matrix<T,ROWS,1> col_vector_at_idx;
        for (size_t row_idx =0; row_idx<ROWS;++row_idx){
            col_vector_at_idx.at(row_idx,0) = at(row_idx,col_idx);
        }             
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::operator*(const Matrix<T,ROWS,COLS>& other)const{
        Matrix<T, ROWS, COLS> result;
        for (size_t i = 0; i < ROWS; ++i) {
            for (size_t j = 0; j < COLS; ++j) {
                result.at(i, j) = at(i,j) * other.at(i, j);
            }
        }
        return result;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS-1,COLS-1> Matrix<T,ROWS,COLS>::getMinor(size_t row_idx, size_t col_idx)const{
        if(row_idx>=ROWS ||col_idx>=COLS){
            throw std::out_of_range(" Index out of range ");
        }

        if( ROWS != COLS){
            throw std::runtime_error(" Matrix is not sqquare");  
        }

        Matrix<T,ROWS-1,COLS-1> minor;
        size_t minor_row = 0;
        for(size_t i=0;i<ROWS;++i){
            if(i==row_idx) continue; // skip to the next iteration of the loop
            size_t minor_col = 0;
            for(size_t j=0;j<COLS;++j){
                if(j==col_idx) continue; 
                minor.at(minor_row,minor_col++)=this->at(i,j);
            }
            ++minor_row;
        }
        return minor;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::identity(){
        static_assert(ROWS==COLS, " Matrix must be square to compute determinant");
        Matrix<T,ROWS,COLS> I;
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                if(i==j){
                    I.at(i,j) = static_cast<T>(1);
                }else{
                    I.at(i,j) = static_cast<T>(0);
                }
            }
        }
        return I;
    }

    template<typename T, size_t ROWS, size_t COLS>
    T Matrix<T,ROWS,COLS>::getDeterminant(){
        static_assert(ROWS==COLS, " Matrix must be square to compute determiant");
        T det = 0;
        if constexpr (ROWS==1){
            det = at(0,0);    
        }else if constexpr (ROWS==2){
            det = at(0,0)*at(1,1)-at(0,1)*at(1,0);
        }else{
            for(size_t col_idx=0;col_idx<COLS;++col_idx){
                T sign = (col_idx%2 == 0) ?1:-1;
                T minor_det = getMinor(0,col_idx).getDeterminant();
                det += sign*at(0,col_idx)*minor_det;
            }
        }
        return det;

    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,COLS,ROWS> Matrix<T,ROWS,COLS>::getTransposed()const{
        Matrix<T,COLS,ROWS> tranpose_matrix;
        for(size_t j=0;j<COLS;++j){
            for(size_t i=0;i<ROWS;++i){
                tranpose_matrix.at(j,i)=this->at(i,j);
            }
        }
        return tranpose_matrix;
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,COLS,ROWS> Matrix<T,ROWS,COLS>::getAdjoint() const{
        Matrix<T,COLS,ROWS> adjoint_matrix;
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                Matrix<T, ROWS-1,COLS-1> minor_matrix = getMinor(i,j);

                T det=minor_matrix.getDeterminant();
                T sign = ((i+j)%2==0)?1:-1;
                T cofactor = sign*minor_matrix.getDeterminant();

                // //transpose the cofactor to obtain the adjugate
                adjoint_matrix.at(j,i) = cofactor;

            }
        }
        return adjoint_matrix;
    }


    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> Matrix<T,ROWS,COLS>::getAdjugate() const{
        return getAdjoint();
    }



    template<typename T, size_t ROWS, size_t COLS>
    template<typename S>
    Matrix<S,ROWS,COLS> Matrix<T,ROWS,COLS>::getInverse(){
        static_assert(ROWS==COLS, " Matrix must be square to compute determinant");
        T determinant = getDeterminant();
        
        assert(static_cast<int>(determinant)!=0);
        
        Matrix<S,ROWS,COLS> inverse;
        
        Matrix<T,ROWS,COLS> Adjoint = getAdjoint();

        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                inverse.at(i,j) =   static_cast<S>(Adjoint.at(i,j)) / static_cast<S>(determinant);
            }
        }

        return inverse;
    }


    /* NON-MEMBER FUNCTIONS */

    //Another matrix with scalar
    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> operator*(T scalar, const Matrix<T,ROWS,COLS>& other) {
        Matrix<T,ROWS,COLS> result;
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                result.at(i,j) = scalar*other.at(i,j);
            }
        }
        return result;
    }  

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> operator/(const Matrix<T,ROWS,COLS>& other, T scalar) {
        Matrix<T,ROWS,COLS> result;
        assert(static_cast<T>(scalar)!=0);
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                result.at(i,j) = other.at(i,j)/scalar;
            }
        }
        return result;
    } 

    // Matrix Multiplication
    template <class T, size_t ROWS1, size_t COLS1, size_t COLS2>
    Matrix<T,ROWS1,COLS2> operator*(const Matrix<T,ROWS1,COLS1>& lhs, const Matrix<T,COLS1,COLS2>& rhs){
        Matrix<T,ROWS1,COLS2> result;
        for(size_t i=0;i<ROWS1;++i){
            for(size_t j=0;j<COLS2;++j){
                result.at(i,j)=0;
                for (size_t k=0;k<COLS1;++k){
                    result.at(i,j) += lhs.at(i,k)*rhs.at(k,j);
                }
            };
        }
        return result;
    }


    template <class S, class T, size_t ROWS, size_t COLS>
    Matrix<S, ROWS, COLS> convert(const Matrix<T, ROWS, COLS>& mat)
    {
        Matrix<S, ROWS, COLS> result;
        for (size_t i = 0; i < ROWS; ++i)
        {
            for (size_t j = 0; j < COLS; ++j)
            {
                result.at(i,j) = mat.at(i,j);
            }
        }
        return result;
    }


    /*As already pointed out in the comments, you have to stop the recursion by defining the specialization Matrix
    Therefore, you must provide a sepcialization template<class T> GenericMatrix<T,1,1> whose GetDeterminant does not call GetMinor
    https://stackoverflow.com/questions/59930894/c-recursive-function-with-template-results-in-zero-sized-array
    */

    template<typename T>
    class Matrix<T,1,1>{
        public:
            Matrix()=default;
            Matrix(std::initializer_list<std::initializer_list<T>> init){
                size_t row = 0;
                for(const auto& row_list: init){
                    //iterate over each element of the row
                    size_t col = 0;
                    for (const auto& value: row_list){
                        at(row,col++)=value;
                    }
                    ++row;
                }
            }

            Matrix<T,1,1>& operator+=(T scalar){
                data[0][0]+=scalar;
                return *this;
            }

            Matrix<T,1,1>& operator/=(T scalar){
                
                data[0][0]/=scalar;
                return *this;
            }

            Matrix<T,1,1>& operator*=(T scalar){    
                data[0][0]*=scalar;
                return *this;
            }
            Matrix<T,1,1> operator+(T scalar) const{
                return Matrix<T,1,1>(*this)+=scalar;
            }

            // Binary arithmetic operators 
            Matrix<T,1,1> operator*(T scalar)const{
                return Matrix<T,1,1>(*this)*=scalar;
            }


            T& at(size_t row, size_t col){
                return data[row][col];
            }
            const T& at(size_t row, size_t col) const{
                return data[row][col];
            }


        private:
            T data[1][1];
    };



}

#endif