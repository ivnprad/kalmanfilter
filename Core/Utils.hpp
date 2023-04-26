#ifndef UTILS_H
#define UTILS_H

#include <ostream>
#include <random>
//#include "onedimension.hpp"
//#include "matrix.hpp"

//forward declarations of matrix.hpp in order to use utils.hpp 
namespace robotics{
    template<typename T,size_t ROWS, size_t COLS>
    class Matrix;

    template<typename T,size_t ROWS, size_t COLS>
    class ConstMatrix;

    template<typename T,size_t ROWS, size_t COLS>
    std::ostream& operator<<(std::ostream& stream, const Matrix<T,ROWS,COLS>& other){
        for(size_t i=0;i<ROWS;++i){
            for(size_t j=0;j<COLS;++j){
                stream << other.at(i,j) << " ";
            }
            stream << "\n";
        }
        return stream;
    }

    template<typename T, size_t ROWS, size_t COLS>
    std::pair<size_t, size_t> getMaxElementIndex(const Matrix<T,ROWS,COLS>& other){
        auto max_it = std::max_element(& other.at(0,0), &other.at(0,0) + ROWS * COLS);
        size_t max_index = std::distance(&other.at(0,0), max_it);
        size_t max_row = max_index / COLS;
        size_t max_col = max_index % COLS;
        return std::make_pair(max_row, max_col);
    }

    template<typename T, size_t ROWS, size_t COLS>
    Matrix<T,ROWS,COLS> power(const Matrix<T,ROWS,COLS>& other, int exponent){
        Matrix<T,ROWS,COLS> result;
        for(size_t row=0; row<ROWS;++row){
            for(size_t col=0;col<COLS;++col){
                result.at(row,col)= pow(other.at(row,col),exponent);
            }
        }
        return result;
    }

    template<typename T, size_t ROWS, size_t COLS>
    constexpr T MatrixTrace(const ConstMatrix<T, ROWS, COLS>& matrix) {
        T trace = 0;
        for (size_t i = 0; i < ROWS; ++i) {
            trace += matrix.atCT(i, i);
        }
        return trace;
    }

    

    template<typename T, size_t ROWS, size_t COLS, size_t row=0, size_t col=0>
    constexpr T MatrixSumUnroll(const ConstMatrix<T,ROWS,COLS>& _matrix){
        if constexpr (row == ROWS) {  // base case: sum of all elements is 0
            return T{};
        }
        else if constexpr (col == COLS - 1) {  // end of row reached, move to next row
            return _matrix.atCT(row, col) + MatrixSumUnroll<T, ROWS, COLS, row + 1, 0>(_matrix);
        }
        else {  // sum current element and move to next column
            return _matrix.atCT(row, col) + MatrixSumUnroll<T, ROWS, COLS, row, col + 1>(_matrix);
        }
    }


}

/// @brief give random value between 0 and 1
/// @return the random value
double myRandom() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return dis(gen);
}

template<typename T=int>
T modulo(T dividend, T divisor) {
    assert(divisor!=0);
    T q = dividend / divisor;
    T p = q * divisor;
    T r = dividend - p;
    if (r < 0) {
        r += divisor;
    }
    return r;
}



/// @brief compile-time square 
/// @tparam T 
/// @param x 
/// @return x*x
template<typename T=double>
constexpr T square(T x) {
    return x * x;
}

constexpr double constexpr_sqrt(double x, double current, double epsilon) {
    return (current - x / current) < epsilon ? current :
           constexpr_sqrt(x, 0.5 * (current + x / current), epsilon);
}

/// @brief square root
/// @param x 
/// @return 
constexpr double sqrt_helper(double x) {
    return x >= 0 && x < 1 ? constexpr_sqrt(x, 1, 0.00001) :
           x >= 1 ? constexpr_sqrt(x, x / 2, 0.00001) :
           std::numeric_limits<double>::quiet_NaN();
}


// pythonysh functions
std::vector<double> range(double x){
    std::vector<double> X;
    for(double i=0;i<x;i++){
        X.push_back(i);
    }
    return X;
}

std::vector<double> randn(size_t nElement){
    std::vector<double> X;
    for(size_t i=0;i<nElement;++i){
        X.push_back(myRandom());
    }
    return X;
}

double mean(const std::vector<double>& other){

    double sum;
    for(const auto& x:other){
        sum+=x;
    }
    return sum/other.size();
}


std::ostream& operator<<(std::ostream& stream, const std::vector<double>& other){
    for(const auto& x: other){
        stream << x << " ";
    }
    stream << "\n";
    return stream;
}

// std::ostream& operator<<(std::ostream& stream, const gaussianS& other){
//     stream << "mean " << other.mean << " variance "<< other.var << "\n";
//     return stream;
// }


std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result.push_back(vec.at(i) * scalar);
    }
    return result;
}

std::vector<double> operator+(const std::vector<double>& vec, double scalar) {
    std::vector<double> result;
    result.reserve(vec.size());
    for (int i = 0; i < vec.size(); i++) {
        result.push_back(vec.at(i) + scalar);
    }
    return result;
}




#endif