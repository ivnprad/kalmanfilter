#ifndef UTILS_H
#define UTILS_H

#include <ostream>
#include <random>
//#include "matrix.hpp"

//forward declarations of matrix.hpp in order to use utils.hpp 
namespace robotics{
    template<typename T,size_t ROWS, size_t COLS>
    class Matrix;

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

}

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




#endif