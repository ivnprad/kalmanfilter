#ifndef TRAIN_H
#define TRAIN_H

#include <random>
#include "matrix.hpp"
#include "Utils.hpp"

using namespace robotics;
using localUnits = int8_t;

template<size_t KERNELCOLS=1>
class Train{
public:
    Train()=delete;
    Train(localUnits track_len, Matrix<double,1,KERNELCOLS> kernel, double sensor_accuracy=0.9 );
    localUnits move(localUnits distance = 1);
    localUnits sense();
    inline localUnits getTrainPos(){return pos; };

private:
    const localUnits m_trackLen;
    const Matrix<double,1,KERNELCOLS> m_kernel;
    const double sensor_accuracy;
    localUnits pos=0;

};


template<size_t KERNELCOLS>
Train<KERNELCOLS>::Train(localUnits track_len, Matrix<double,1,KERNELCOLS> kernel, double sensor_accuracy)
  : m_trackLen(track_len), m_kernel(kernel), sensor_accuracy(sensor_accuracy) {}


/// @brief .
//////////""" move in the specified direction with some small chance of error"""
///////// the idea is to displace the pos to any of the positions meant the by the kernel 
//////// e.g. if the kernel is [0.1,0.8,0.1] and myRandom gives 0.9 then it witha distance of 1 
/////// it should no alter the pos of the sensor
/////// if myRandom gives 0.2 it should displace -1 + distance 
//////  if my Randome gives 1.0 it should displace +1 + distance 
//////////add noise here if the sensor is perfect that means m_kernel = [1.] and myRandom returns maximum 1
////////// offset should be zero. Other offset should increment until the accumulate of probabilites of m_kenel
////////// is higher than rnd 
////////// The IDEA 
/// @tparam KERNELCOLS 
/// @param distance 
/// @return 
template<size_t KERNELCOLS>
localUnits Train<KERNELCOLS>::move(localUnits distance){
    pos+=distance;
    const double rnd = myRandom();
    double s=0.0;
    localUnits offset = -(static_cast<localUnits>(KERNELCOLS)-1)/2;
    for(size_t idx = 0; idx<KERNELCOLS;++idx){
        s+=m_kernel.at(0,idx);
        if (s>=rnd){
            break;
        }
        offset += 1;
    }
    pos = modulo(static_cast<localUnits>(pos+offset),static_cast<localUnits>(m_trackLen));
    return pos;

}

template<size_t KERNELCOLS>
localUnits Train<KERNELCOLS>::sense(){
    double itsPos = pos;
    if (myRandom()>sensor_accuracy){
        if(myRandom()>0.5){
            itsPos+=1;
        }else{
            itsPos-=1;
        }
    }
    return itsPos;
}


#endif 