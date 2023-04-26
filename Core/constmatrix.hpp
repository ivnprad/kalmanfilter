#ifndef CONSTMATRIX_HPP
#define CONSTMATRIX_HPP

#include <iostream>
#include <cmath>
#include <initializer_list>

namespace robotics{

    template<typename T, size_t ROWS, size_t COLS>
    class ConstMatrix {
    public:
        // Default constructor
        ConstMatrix() = default;

        // Constructor that accepts an initializer list
        constexpr ConstMatrix(std::initializer_list<std::initializer_list<T>> list) {
            size_t rowIdx = 0;
            for (const auto& row : list) {
                size_t colIdx = 0;
                for (const auto& val : row) {
                    _data[rowIdx][colIdx++] = val;
                }
                ++rowIdx;
            }
        }

        // Member function to access elements
        T& at(size_t row, size_t col) {
            return _data[row][col];
        }

        // Member function to access elements at compile-time
        constexpr T atCT(size_t row, size_t col) const {
            return _data[row][col];
        }

    private:
        T _data[ROWS][COLS];
    };

}

#endif 