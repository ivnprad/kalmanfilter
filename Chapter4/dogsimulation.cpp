#include "dogsimulation.hpp"
//#include "Utils.hpp"
#include <cmath>
#include <random>

Dogsimulation::Dogsimulation(double x0,double _velocity,double measurement_var,double process_var):
    x(x0),velocity(_velocity){
        meas_std = sqrt(measurement_var);
        process_std=sqrt(process_var);
        randomNumber = []() { 
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            return dis(gen);
        };
}

void Dogsimulation::move(double dt){
    const double dx = velocity + randomNumber()*process_std;
    x+=dx*dt; 
} 

double Dogsimulation::sense_position(){
    return x + randomNumber()*meas_std;
}

double Dogsimulation::move_and_sense(){
    move();
    return  sense_position(); 
}