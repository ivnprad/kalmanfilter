#include "possensor.hpp"

#include <cmath>
#include <random>
//#include "../Core/Utils.hpp"



PosSensor::PosSensor(std::vector<double> _pos, std::vector<double> _vel,double _noise_std):
    pos(_pos),vel(_vel), noise_std(_noise_std){
        
    randomNumber = [](const double lowerlimit=-1.0, const double upperlimit=1.0) { 
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(lowerlimit, upperlimit); // to get negative numbers -1.0 
            return dis(gen);
        };
    
}

DynamicMatrix<double> PosSensor::read(){
    pos.at(0)+=vel.at(0);
    pos.at(1)+=vel.at(1);

    return DynamicMatrix<double>({{pos.at(0)+randomNumber()*noise_std},
                                    {pos.at(1)+randomNumber()*noise_std}});


}