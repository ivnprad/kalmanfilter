#ifndef ONEDIMENSION_HPP
#define ONEDIMENSION_HPP

#include <vector>
#include <cmath>
//#include "probabilities.hpp"
#include "../Core/matplotlibcpp.h"
#include "dogsimulation.hpp"

namespace plt = matplotlibcpp;

    // CHAPTER 4

//void plot_gaussian_pdf(double mean, double variance,std::vector<double> xlim, std::vector<double> ylim){

double gaussianDensity(double x, double mu, double sigma) {
    double coefficient = 1.0 / (sigma * sqrt(2 * M_PI));
    double exponent = -0.5 * pow((x - mu) / sigma, 2);
    return coefficient * exp(exponent);
}

void plot_gaussian_pdf(double mean, double variance,std::vector<double> xlim){
    double std = sqrt(variance);

    std::vector<double> xs, ys;
    double step = (xlim.at(1) - xlim.at(0)) / 1000.0;
    for (double i = xlim.at(0); i < xlim.at(1); i += step) {
        xs.push_back(i);
        const double pdf_at_i = gaussianDensity(i,mean,variance);
        ys.push_back(pdf_at_i);
    }

    // Create a new plot
    plt::plot(xs, ys,"k-");

    //plt::ylim(0,1.0);
    plt::xlim(6, 14);
    plt::ylim(0,1);
    // Show the plot
    //plt::show();
}

struct gaussianS {
        double mean;
        double var;
        gaussianS(double _mean, double _variance) : mean(_mean),var(_variance){}
};

gaussianS operator*(const gaussianS& lhs,const gaussianS& rhs){
    const double mean = (lhs.mean*rhs.var+rhs.mean*lhs.var)/(lhs.var+rhs.var);
    const double var = rhs.var*lhs.var/(lhs.var+rhs.var);
    return gaussianS(mean,var);
}

void plot_products(gaussianS lhs, gaussianS rhs, std::vector<double> xlim){
    plot_gaussian_pdf(lhs.mean, lhs.var, xlim);
    plot_gaussian_pdf(rhs.mean, rhs.var, xlim);
    const auto product = lhs*rhs;
    plot_gaussian_pdf(product.mean,product.var,xlim );
    plt::show();
}

gaussianS operator+(const gaussianS& lhs,const gaussianS& rhs){
    return gaussianS(lhs.mean+rhs.mean,lhs.var+rhs.var);
}



std::ostream& operator<<(std::ostream& stream, const gaussianS& other){
    stream << "mean " << other.mean << " variance " << other.var << "\n";
    return stream;
}

gaussianS predict(gaussianS pos, gaussianS movement){
    return pos+movement;
}

gaussianS updateODKF(gaussianS likelihood, gaussianS prior){
    return likelihood*prior;
}

#endif 