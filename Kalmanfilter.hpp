#ifndef SUPERVISEDLEARNING_H
#define SUPERVISEDLEARNING_H


#include "matrix.hpp"
#include <random>


/**
 * ROWS  Training Examples
 * COLS  Features
**/

using namespace robotics;

/******  FEATURE SCALING ************
 * to choose learning rate start with a small number and keep increasing 3X  
 */
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


/*weight GOLD */
template<typename T, size_t ROWS>
Matrix<T,ROWS,1> g_h_filter(const Matrix<T,ROWS,1>& other, const T& X0, T dx, T g, T h, T dt){
    /*
    """
    Performs g-h filter on 1 state variable with a fixed g and h.
    ’data’ contains the data to be filtered.
    ’x0’ is the initial value for our state variable
    ’dx’ is the initial change rate for our state variable 
    ’g’ is the g-h’s g scale factor
    ’h’ is the g-h’s h scale factor
    ’dt’ is the length of the time step 
    """
     */
    size_t rows = other.getRows();
    Matrix<T,ROWS,1> estimates;
    T x = X0;

    for(size_t row =0; row<rows;++row){
        //1 .- measurement
        const T& z = other.at(row,0);

        //2.- predic new position using systems dynamic model  STATE EXTRAPOLATION EQUATIONs
        const T x_estimate = x + dx*dt; //Xn+1,n
        dx = dx; // 

        //3.- update STATES  STATE UPDATE EQUATION
        const T innovation = z-x_estimate; // innovation or residual Zn
        dx = dx + h*(innovation/dt); // UPDATE DX
        x = x_estimate + g*innovation; // Xn,n-1 ( this comes from step 2 with unit delay ), Zn , g -> Xn,n

        estimates.at(row,0) = x;

    } 
    return estimates;
}

template<typename T, size_t ROWS> // instead of count 
Matrix<T,ROWS,1> gen_data(const T& X0, T dx,T noisy_factor){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    Matrix<double,ROWS,1> generatedData;

    for(size_t i=0;i<ROWS;++i){
        const double randomValue = distribution(gen);
        generatedData.at(i,0) = X0+dx*i+randomValue*noisy_factor;
    }
    
    return generatedData;

}


#endif 
