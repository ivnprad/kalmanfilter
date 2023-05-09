#ifndef MULTIVARIATEGAUSSIANS_HPP
#define MULTIVARIATEGAUSSIANS_HPP

//#include "../Core/matrix.hpp"
#include "../Core/dynamicMatrix.hpp"
#include <vector>
#include <algorithm>

//using namespace robotics;

template<typename T=double>
DynamicMatrix<T> covariance(const std::vector<T>& X, const std::vector<T>& Y, const uint8_t bias=0){
    assert(Y.size()==X.size());
    T meanX=0;
    T meanY=0;
    T varianceX=0;
    T varianceY=0;
    T covarianceXY=0;

    
    const T nELE = static_cast<T>(X.size());

    for(size_t idx=0;idx<nELE;++idx){
        meanX+=X[idx];
        meanY+=Y[idx];
    }
    meanX/=nELE;
    meanY/=nELE;

    for(size_t idx=0;idx<nELE;++idx){
        varianceX+=pow(X[idx]-meanX,2);
        varianceY+=pow(Y[idx]-meanY,2);
        covarianceXY+=(X[idx]-meanX)*(Y[idx]-meanY);
    }
    // biased estimator
    const T correction = bias==0?1:0;
    varianceX/=nELE- correction;
    varianceY/=nELE-correction;
    covarianceXY/=nELE-correction;

    DynamicMatrix cov{{varianceX,covarianceXY},{covarianceXY,varianceY}};

    return cov;
}


template<typename T=double>
T mutlivariateGaussian(const std::vector<T>& X,const std::vector<T>& mu, DynamicMatrix<T> covariance) {
    const size_t nVARIABLES=2;
    assert(covariance.rows()==nVARIABLES);
    assert(X.size()==nVARIABLES);
    assert(covariance.rows()==covariance.cols());
    assert(X.size()==mu.size());
    
    const DynamicMatrix<T> x_mu({{X[0]-mu[0]},{X[1]-mu[1]}});
    const DynamicMatrix<T> x_mu_T = x_mu.getTransposed();
    const T determinant = covariance.getDeterminant();
    const DynamicMatrix<T> itsInverse = covariance.getInverse();
    
    const double sigma = sqrt(determinant);
    const double coefficient = 1.0/(sigma*2*M_PI); // sqrt(pow(2*M_PI,2))
    //std::cout << " coefficient "<< coefficient<< std::endl;
    const double exponent =  0.5*(((x_mu_T^itsInverse)^x_mu)[0][0]);
    //std::cout << " exponent "<< exponent << std::endl;
    return coefficient*exp(-exponent);

}

template<typename T=double>
std::pair<DynamicMatrix<T>,std::vector<T>>  calculate_eigenvectors(DynamicMatrix<T> covariance_matrix) {
    
    DynamicMatrix<T> eigenvectors({{0.0,0.0},{0.0,0.0}});

    // Calculate the eigenvalues
    double trace = covariance_matrix(0,0) + covariance_matrix(1,1);
    double determinant = covariance_matrix.getDeterminant();
    double discriminant = sqrt(trace * trace - 4 * determinant);
    double lambda1 = (trace + discriminant) / 2;
    double lambda2 = (trace - discriminant) / 2;

    std::vector<T> eigenvalues;
    eigenvalues.emplace_back(lambda1);
    eigenvalues.emplace_back(lambda2);

    // Calculate the eigenvectors
    if (covariance_matrix(0,1) == 0) {
        eigenvectors(0,0) = 1.0;
        eigenvectors(1,1) = 1.0;
    }
    else {
        eigenvectors(0,0) = -(covariance_matrix(1,1)-lambda1) ;
        eigenvectors(1,0) = covariance_matrix(0,1);
        eigenvectors(0,1) = -covariance_matrix(0,1);
        eigenvectors(1,1) = (covariance_matrix(0,0)-lambda2);
    }
    // Normalize the eigenvectors
    double length1 = sqrt(eigenvectors(0,0) * eigenvectors(0,0) + eigenvectors(1,0) * eigenvectors(1,0));
    double length2 = sqrt(eigenvectors(0,1) * eigenvectors(0,1) + eigenvectors(1,1) * eigenvectors(1,1));
    eigenvectors(0,0) /= length1;
    eigenvectors(1,0) /= length1;
    eigenvectors(0,1) /= length2;
    eigenvectors(1,1) /= length2;

    return std::make_pair(eigenvectors,eigenvalues);
}


void plot_ellipse(double lambda1, double lambda2, double v1_x, double v1_y, double v2_x, double v2_y, double stdDevMultiples =1)
{
    double theta = atan2(v1_y, v1_x) + M_PI / 2.0;
    std::cout << " angle " << theta*180/M_PI << std::endl;
    //double angle = atan2(v2_y, v2_x);

    double Cx = 0;
    double Cy = 0;
    double Rminor = 0; 
    double Rmajor = 0; 
    if( lambda1>lambda2){
        Rminor = sqrt(lambda2)*stdDevMultiples;
        Rmajor = sqrt(lambda1)*stdDevMultiples; 
    }else{
        Rminor = sqrt(lambda1)*stdDevMultiples; 
        Rmajor = sqrt(lambda2)*stdDevMultiples; 
    }

    const size_t nELE = 1000;
    std::vector<double> x(nELE+1), y(nELE+1);
    for (int i = 0; i <= nELE; ++i)
    {
        double alpha = i * 2*M_PI / static_cast<double>(nELE);
        x.at(i) = Cx + Rmajor*cos(alpha)*cos(theta)-Rminor*sin(alpha)*sin(theta);
        y.at(i) = Cy + Rmajor*cos(alpha)*sin(theta)+Rminor*sin(alpha)*cos(theta);
    }    

    plt::plot(x, y);
    plt::set_aspect_equal();
    plt::show();
}



#endif 