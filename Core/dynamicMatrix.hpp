#ifndef DYNAMICMATRIX_HPP
#define DYNAMICMATRIX_HPP

#include <vector>
#include <iostream>
#include <type_traits>
#include "../Core/Utils.hpp"

using namespace std;

//template <typename T>
template <typename T=double>
class DynamicMatrix {
public:
    // Constructor that creates a matrix with n rows and m columns
    DynamicMatrix(int rows, int cols):data(rows, std::vector<T>(cols)){}

    // Constructor that creates a matrix with a given initializer list
    DynamicMatrix(std::initializer_list<std::initializer_list<T>> init);

    // Accessors
    T& operator()(int row, int col);
    const T& operator()(int row, int col) const ;
    inline int rows() const { return data.size();}
    inline int cols() const { return data[0].size();}
    std::vector<T>& operator[](int row);
    const std::vector<T>& operator[](int row) const;
    std::vector<T>& getRow(int row);    
    const std::vector<T>& getRow(int row) const;
    std::vector<T> getCol(int col);    
    const std::vector<T> getCol(int col) const;
    
    //comparators
    bool operator==(const DynamicMatrix<T>& other) const;
    bool operator!=(const DynamicMatrix<T>& other) const;

    //compound assigment operators
    DynamicMatrix<T>& operator+=(const DynamicMatrix<T>& other);
    DynamicMatrix<T>& operator-=(const DynamicMatrix<T>& other);
    DynamicMatrix<T>& operator/=(const DynamicMatrix<T>& other);
    DynamicMatrix<T>& operator*=(const DynamicMatrix<T>& other); // element wise multiplication 
    DynamicMatrix<T>& operator*=(T scalar);
    DynamicMatrix<T>& operator/=(T scalar);
    DynamicMatrix<T>& operator+=(T scalar);
    DynamicMatrix<T>& operator-=(T scalar);

    //binary arithmetic operators
    DynamicMatrix<T> operator*(T scalar) const;
    DynamicMatrix<T> operator+(T scalar) const;
    DynamicMatrix<T> operator-(T scalar) const;
    DynamicMatrix<T> operator+(const DynamicMatrix<T>& other)const;
    DynamicMatrix<T> operator-(const DynamicMatrix<T>& other)const;
    DynamicMatrix<T> operator/(const DynamicMatrix<T>& other)const; //#TODO this might be wrong implemented
    DynamicMatrix<T> operator^(const DynamicMatrix<T>& other)const; //MATRIX MULTIPLICATION
    DynamicMatrix<T> operator*(const DynamicMatrix<T>& other)const; // element wise matrix multiplication 

    // Linear Algebra
    DynamicMatrix<T> getMinor(size_t row_idx, size_t col_idx)const;
    DynamicMatrix<T> identity();
    //static const Matrix<T,ROWS,COLS> GetRotMatrix(const Quaternion& q_rot); // later
    T getDeterminant();
    DynamicMatrix<T> getTransposed() const;
    DynamicMatrix<T> getAdjoint() const;
    DynamicMatrix<T> getAdjugate() const;

    template <typename S = double>
    DynamicMatrix<S> getInverse();    

private:
    std::vector<std::vector<T>> data;
};

template<typename T>
DynamicMatrix<T>::DynamicMatrix(std::initializer_list<std::initializer_list<T>> init) {
    for (const auto& row : init) {
        data.emplace_back(row);
    }
}


template<typename T>
T& DynamicMatrix<T>::operator()(int row, int col) {
    if (row < 0 || row >= data.size() || col < 0 || col >= data[0].size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row][col];
}

template<typename T>
const T& DynamicMatrix<T>::operator()(int row, int col) const {
    if (row < 0 || row >= data.size() || col < 0 || col >= data[0].size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row][col];
}



template<typename T>
std::vector<T>& DynamicMatrix<T>::operator[](int row) {
    if (row < 0 || row >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row];
}

template<typename T>
const std::vector<T>& DynamicMatrix<T>::operator[](int row) const {
    if (row < 0 || row >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row];
}

template<typename T>
std::vector<T>& DynamicMatrix<T>::getRow(int row){
    if (row < 0 || row >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row];
}

template<typename T>
const std::vector<T>& DynamicMatrix<T>::getRow(int row) const{
    if (row < 0 || row >= data.size()) {
        throw std::out_of_range("Index out of range");
    }
    return data[row];
}

template<typename T>
std::vector<T> DynamicMatrix<T>::getCol(int col){
    const size_t ROWS=this->data.size();
    std::vector<T> column(ROWS);
    for(size_t row =0; row<ROWS; ++row){
        column.at(row)=data[row][col];
    }
    return column;
}

template<typename T>    
const std::vector<T> DynamicMatrix<T>::getCol(int col) const{
    const size_t ROWS=rows();
    std::vector<T> column(ROWS);
    for(size_t row =0; row<ROWS; ++row){
        column.at(row)=data[row][col];
    }
    return column;    
}



//bool operator
template<typename T>
bool DynamicMatrix<T>::operator==(const DynamicMatrix<T>& other) const{
    if(this==&other){// comparing with itself
        return true;
    }
    const size_t ROWS=rows();
    const size_t COLS=cols();
    return std::equal(&data[0][0],&data[ROWS-1][COLS-1]+1,&other.data[0][0]);
}

template<typename T>
bool DynamicMatrix<T>::operator!=(const DynamicMatrix<T>& other) const{
    return !(*this==other);
}

//compound assigment operators
template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator+=(const DynamicMatrix<T>& other){
    const size_t ROWS=rows();
    const size_t COLS=cols();    
    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            (row,col) += other(row,col);
        }
    }
    return *this;    
}


template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator-=(const DynamicMatrix<T>& other){
    const size_t ROWS=rows();
    const size_t COLS=cols();    
    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            (row,col) -= other(row,col);
        }
    }
    return *this;       
}

template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator/=(const DynamicMatrix<T>& other){
    const size_t ROWS=rows();
    const size_t COLS=cols();        
    for (size_t row = 0; row < ROWS; ++row) {
        for (size_t col = 0; col < COLS; ++col) {
            (row,col) /= other(row,col);
        }
    }
    return *this;    
}

template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator*=(const DynamicMatrix<T>& other){
    const size_t ROWS=rows();
    const size_t COLS=cols();        
    for (size_t row = 0; row < ROWS; ++row) {
        for (size_t col = 0; col < COLS; ++col) {
            (row,col) *= other(row,col);
        }
    }
    return *this;    
}


template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator*=(T scalar){
    const size_t ROWS=rows();
    const size_t COLS=cols();     
    std::for_each(&(data[0][0]),&(data[ROWS-1][COLS-1])+1,[=](auto& element){
        element *= scalar;
    });
    return *this;    
}

template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator/=(T scalar){
    assert(static_cast<int>(scalar)!=0);
    T new_scalar = 1/scalar;
    (*this)*=new_scalar;
    return *this;
}

template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator+=(T scalar){
    const size_t ROWS=rows();
    const size_t COLS=cols();    
    std::for_each(&(data[0][0]),&(data[ROWS-1][COLS-1])+1,[=](auto& element){
        element += scalar;
    });
    return *this;
}

template<typename T>
DynamicMatrix<T>& DynamicMatrix<T>::operator-=(T scalar){
    const size_t ROWS=rows();
    const size_t COLS=cols();    
    std::for_each(&(this->data[0][0]),&(this->data[ROWS-1][COLS-1])+1,[=](auto& element){
        element -= scalar;
    });
    return *this;
}

// Binary arithmetic operators 
template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator*(T scalar)const{
    return DynamicMatrix<T>(*this)*=scalar;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator+(T scalar) const{
    return DynamicMatrix<T>(*this)+=scalar;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator-(T scalar) const{
    return DynamicMatrix<T>(*this)-=scalar;
}    

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator+(const DynamicMatrix<T>& other) const{
    return DynamicMatrix<T>(*this)+=other;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator-(const DynamicMatrix<T>& other) const{
    return DynamicMatrix<T>(*this)-=other;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator/(const DynamicMatrix<T>& other)const{
    return DynamicMatrix<T>(*this)/=other;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::operator^(const DynamicMatrix<T>& other)const{
    const size_t M=rows();
    const size_t N = cols();
    const size_t rhsN = other.rows();
    const size_t P=other.cols();

    assert(N==rhsN);

    DynamicMatrix<T> result(M,P);
    for(size_t row=0; row<M; ++row){
        for(size_t col=0; col<P; ++col){
            result(row,col) = getRow(row) * other.getCol(col);
        }

    }

    return result;
}


template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::getMinor(size_t row_idx, size_t col_idx)const{
    const size_t ROWS=rows();
    const size_t COLS=cols();

    // if(row_idx>=ROWS ||col_idx>=COLS){
    //     throw std::out_of_range(" Index out of range ");
    // }

    if( ROWS != COLS){
        throw std::runtime_error(" Matrix is not sqquare");  
    }

    DynamicMatrix<T> minor(ROWS-1,COLS-1);
    size_t minor_row = 0;
    for(size_t i=0;i<ROWS;++i){
        if(i==row_idx) continue; // skip to the next iteration of the loop
        size_t minor_col = 0;
        for(size_t j=0;j<COLS;++j){
            if(j==col_idx) continue;
            minor(minor_row,minor_col++)=data[i][j];
        }
        ++minor_row;
    }
    return minor;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::identity(){
    const size_t ROWS=data.size();
    const size_t COLS=data[0].size();
    assert(ROWS==COLS);
    DynamicMatrix<T> I;
    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            if(row==col){
                I(row,col) = static_cast<T>(1);
            }else{
                I(row,col) = static_cast<T>(0);
            }
        }
    }
    return I;
}

template<typename T>
T DynamicMatrix<T>::getDeterminant(){
    const size_t ROWS=rows();
    const size_t COLS=cols();
    assert(ROWS==COLS);
    T det = 0;
    if (ROWS==1){
        det = data[0][0];    
    }else if (ROWS==2){
        det = data[0][0]*data[1][1]-data[0][1]*data[1][0];
    }else{
        for(size_t col_idx=0;col_idx<COLS;++col_idx){
            T sign = (col_idx%2 == 0) ?1:-1;
            T minor_det = getMinor(0,col_idx).getDeterminant();
            det += sign*data[0][col_idx]*minor_det;
        }
    }
    return det;

}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::getTransposed()const{
    const size_t ROWS=rows();
    const size_t COLS=cols();

    DynamicMatrix<T> tranpose_matrix(COLS,ROWS);
    for(size_t col=0;col<COLS;++col){ 
        for(size_t row=0;row<ROWS;++row){
            tranpose_matrix(col,row)=data[row][col];
        }
    }
    return tranpose_matrix;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::getAdjoint() const{
    const size_t ROWS=rows();
    const size_t COLS=cols();    
    DynamicMatrix<T> adjoint_matrix(COLS,ROWS);
    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            DynamicMatrix<T> minor_matrix = getMinor(row,col);

            //T det=minor_matrix.getDeterminant();
            T sign = ((row+col)%2==0)?1:-1;
            T cofactor = sign*minor_matrix.getDeterminant();

            // //transpose the cofactor to obtain the adjugate
            adjoint_matrix(col,row) = cofactor;

        }
    }
    return adjoint_matrix;
}

template<typename T>
DynamicMatrix<T> DynamicMatrix<T>::getAdjugate() const{
    return getAdjoint();
}

template<typename T>
template<typename S>
DynamicMatrix<S> DynamicMatrix<T>::getInverse(){
    const size_t ROWS=rows();
    const size_t COLS=cols();        
    assert(ROWS==COLS);
    T determinant = getDeterminant();
    
    assert(static_cast<int>(determinant)!=0);
    
    DynamicMatrix<S> inverse(ROWS,COLS);
    DynamicMatrix<T> Adjoint = getAdjoint();

    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            inverse(row,col) =   static_cast<S>(Adjoint(row,col)) / static_cast<S>(determinant);
        }
    }

    return inverse;
}


// // Output stream operator for Matrix
template <typename T>
std::ostream& operator<<(std::ostream& os, const DynamicMatrix<T>& matrix) {
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            os << static_cast<T>(matrix(i, j)) << " ";
        }
        os << std::endl;
    }
    return os;
}



#endif