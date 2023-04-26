#ifndef PROBABILITIES_HPP
#define PROBABILITIES_HPP

//#include <array>
#include <cmath>
#include "../Core/matrix.hpp"
#include "../Core/constmatrix.hpp"
#include "../Core/Utils.hpp"

using namespace robotics;

    // CHAPTER 3

template<size_t COLS, size_t ROWS =1, typename T = double>
Matrix<T,ROWS,COLS> generateRandomNumbers(){
    Matrix<T,ROWS,COLS> matrixOfRandomNumbers(0);
    for(size_t row=0;row<ROWS;++row){
        for(size_t col = 0; col<COLS;++col){
            matrixOfRandomNumbers.at(row,col) =  myRandom();
        }
    }
    return matrixOfRandomNumbers;
}

template<typename T, size_t ROWS, size_t COLS>
constexpr T mean(const ConstMatrix<T,ROWS,COLS>& _matrix){
    constexpr T totalNumberOfElements = ROWS * COLS;
    //return MatrixSumUnroll<T,ROWS,COLS>(_matrix)/totalNumberOfElements;
    return MatrixSumUnroll(_matrix)/totalNumberOfElements;
}


template<typename T, size_t ROWS, size_t COLS, size_t row=0, size_t col=0>
constexpr T square(const ConstMatrix<T,ROWS,COLS>& _matrix,const T value =0){
    
    if constexpr (row == ROWS) {  // base case: sum of all elements is 0
        return T{};
    }
    else if constexpr (col == COLS - 1) {  // end of row reached, move to next row
        return square(_matrix.atCT(row, col)-value) + square<T, ROWS, COLS, row + 1, 0>(_matrix,value);
    }
    else {  // sum current element and move to next column
        return square(_matrix.atCT(row, col)-value) + square<T, ROWS, COLS, row, col + 1>(_matrix,value);
    }
}


template<typename T, size_t ROWS, size_t COLS>
constexpr T variance(const ConstMatrix<T,ROWS,COLS>& _matrix){

    T itsMean = mean(_matrix);
    T itsSquare = square(_matrix,itsMean);
    constexpr T totalNumberOfElements = ROWS * COLS;

    return itsSquare/totalNumberOfElements;
}

template<typename T, size_t ROWS, size_t COLS>
constexpr T standardDeviation(const ConstMatrix<T,ROWS,COLS>& _matrix){

    T stdDev = variance(_matrix);

    return sqrt_helper(stdDev);
}


template<typename T, size_t ROWS, size_t COLS>
Matrix<T,1,COLS> meanColumnwise(const Matrix<T,ROWS,COLS>& other) {
    Matrix<T,1,COLS> itsMean;
    for(size_t j=0;j<COLS;++j){
        T sum=0;
        for(size_t i=0;i<ROWS;++i){
                sum+= other.at(i,j);
        }
        itsMean.at(0,j) = sum/static_cast<T>(ROWS);
    }
    return itsMean;
}  

template<typename T, size_t ROWS, size_t COLS>
Matrix<T,1,COLS> standardDeviationColumnwise(const Matrix<T,ROWS,COLS>& other){

    //, const Matrix<T,1,COLS>& meanColumnwise_) {
    Matrix<T,1,COLS> meanColumnwise_ = meanColumnwise(other);
    Matrix<T,1,COLS> itsStandardDeviation;
    for(size_t j=0;j<COLS;++j){
        T sum_of_squares=0;
        for(size_t i=0;i<ROWS;++i){
            sum_of_squares +=  std::pow(other.at(i,j)-meanColumnwise_.at(0,j),2);
        }
        itsStandardDeviation.at(0,j) =  std::sqrt(sum_of_squares/static_cast<T>(ROWS-1));
        //The factor of N−1 is called Bessel's correction -> ROWS-1

    }

    return itsStandardDeviation;
}

template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> gaussian(Matrix<T,ROWS,COLS> X, double mu, double variance) {
    Matrix<T,ROWS,COLS> itsGaussian;
    const double sigma = sqrt(variance);
    const double coefficient = 1.0 / (sigma * sqrt(2 * M_PI));

    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            const double exponent = -0.5 * pow((X.at(row,col) - mu) / sigma, 2);
            itsGaussian.at(row,col)=coefficient * exp(exponent);
        }
    }

    return itsGaussian;
}

template<typename T, size_t COLS>
std::pair<T,T> mean_var(const Matrix<T,1,COLS>& P){
    size_t counter=0;
    Matrix<T,1,COLS> p = P;
    Matrix<T,1,COLS> x(0);
    p = p/p.sum();

    for(size_t col=0;col<COLS;++col){
        x.at(0,col) = col;
    }
    
    // std::cout << " p " << p << std::endl;
    // std::cout << " x " << x << std::endl;

    const T mean = (p*x).sum();
    const T variance = (power(x-mean,2)*p).sum();


    return std::pair<T,T>(mean,variance);
}

// double gaussianDensity(double x, double mu, double sigma) {
//     double coefficient = 1.0 / (sigma * sqrt(2 * M_PI));
//     double exponent = -0.5 * pow((x - mu) / sigma, 2);
//     return coefficient * exp(exponent);
// }

double gaussianIntegral(double lowerLimit, double upperLimit, double mean, double standarDeviation) {
    double a_normalized = (lowerLimit - mean) / (standarDeviation * sqrt(2.0));
    double b_normalized = (upperLimit - mean) / (standarDeviation * sqrt(2.0));
    double integral = 0.5 * (erf(b_normalized) - erf(a_normalized));
    return integral;
}


template<typename T, size_t COLS>
void pollAnalysisAge(Matrix<T,1,COLS> poll){
    const std::pair<T,T> _MeanVar = mean_var(poll);
    const T _Mean = _MeanVar.first;
    const T _Variance = _MeanVar.second;
    const T _StdDev = sqrt(_Variance);
    std::cout << "mean " << _Mean << " variance " << _Variance << std::endl;
    std::cout << " standardeviation " << _StdDev << std::endl;

    std::cout<< " assuming 1-10 ->1-100 years the mean age interested is " 
    <<  _Mean*10<< " years old" << std::endl;

    const T lowerLimit =  _Mean-1*_StdDev;
    const T upperLimit = _Mean+1*_StdDev;

    std::cout<< " 65% of the people is within " << lowerLimit*10 << " - " << upperLimit*10 << " " << std::endl;

}






#endif 