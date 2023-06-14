#ifndef POSSENSOR_HPP
#define POSSENSOR_HPP

#include <vector>
#include "../Core/dynamicMatrix.hpp"

class PosSensor{
    public:
        PosSensor(std::vector<double> _pos = std::vector<double>(0,0), std::vector<double> _vel = std::vector<double>(0,0),double noise_std=1.0);
        DynamicMatrix<double> read();
    private:
        std::vector<double> vel;
        double noise_std;
        std::vector<double> pos;
        std::function<double(void)> randomNumber;
};

#endif 
