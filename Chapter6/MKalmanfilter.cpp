#include "MKalmanFilter.hpp"

#include <cmath>
#include <limits>

MKalmanFilter::MKalmanFilter( int _dim_x, int _dim_z, int _dim_u)
    :dim_x(_dim_x),
    dim_z(_dim_z),
    dim_u(_dim_u),
    _I(DynamicMatrix<double>::IdentityMatrixCreator<>::create(dim_x))
    {
    x = DynamicMatrix(dim_x,1);
    P = DynamicMatrix<double>::IdentityMatrixCreator<>::create(dim_x); // uncertainty covariance
    Q = DynamicMatrix<double>::IdentityMatrixCreator<>::create(dim_x); // process uncertainty
    F = DynamicMatrix<double>::IdentityMatrixCreator<>::create(dim_x); // state transition matrix 
    R = DynamicMatrix<double>::IdentityMatrixCreator<>::create(dim_z);
    H = DynamicMatrix(dim_z,dim_x);
    z = DynamicMatrix(dim_z,1);

    //# these will always be a copy of x,P after predict() is called


    K = DynamicMatrix(dim_x,dim_z);
    y = DynamicMatrix(dim_z,1);
    S = DynamicMatrix(dim_z,dim_z);
    SI = DynamicMatrix(dim_z,dim_z);

    //
    x_prior = x;
    P_prior = P;
    x_post = x;
    P_post = P;
    _alpha_sq = 1.0;


    log_likelihood = log(std::numeric_limits<double>::min());
    likelihood = std::numeric_limits<double>::min();
    _mahalanobis = 0.0;  // Initialized to 0, similar to None in Python



    
}

int MKalmanFilter::getDimX() const{
    return dim_x;
}

int MKalmanFilter::getDimZ() const{
    return dim_z;
}

DynamicMatrix<double> MKalmanFilter::getX()const{
    return x;
}

DynamicMatrix<double> MKalmanFilter::getP()const{
    return P;
}

DynamicMatrix<double> MKalmanFilter::getX_PRIOR() const{
    return x_prior;
}

DynamicMatrix<double> MKalmanFilter::getP_PRIOR() const{
    return P_prior;
}

DynamicMatrix<double> MKalmanFilter::getX_POST() const{
    return x_post;
}

DynamicMatrix<double> MKalmanFilter::getP_POST() const{
    return P_post;
}

DynamicMatrix<double> MKalmanFilter::getQ()const{
    return Q;
}

DynamicMatrix<double> MKalmanFilter::getR()const{
    return R;
}

DynamicMatrix<double> MKalmanFilter::getF()const{
    return F;
}

DynamicMatrix<double> MKalmanFilter::getH()const{
    return H;
}

DynamicMatrix<double> MKalmanFilter::getK() const{
    return K;
}

DynamicMatrix<double> MKalmanFilter::getY() const{
    return y;
}

DynamicMatrix<double> MKalmanFilter::getS() const{
    return S;
}

DynamicMatrix<double> MKalmanFilter::getSI() const{
    return SI;
}


DynamicMatrix<double> MKalmanFilter::getM() const{
    return M;
}

DynamicMatrix<double> MKalmanFilter::getB() const{
    return B;
}

DynamicMatrix<double> MKalmanFilter::getZ() const{
    return z;
}


double MKalmanFilter::getLogLikeihood() const{
    return log_likelihood;
}

double MKalmanFilter::getLikeihood() const{
    return likelihood;
}

double MKalmanFilter::getMahalanobis() const{
    return _mahalanobis;
}

double MKalmanFilter::getAlpha() const{
    return _alpha_sq;
}

void MKalmanFilter::setX(const DynamicMatrix<double>& other){
    F=other;
}
void MKalmanFilter::setP(const DynamicMatrix<double>& other){
    P=other;
}
void MKalmanFilter::setQ(const DynamicMatrix<double>& other){
    Q=other;
}
void MKalmanFilter::setH(const DynamicMatrix<double>& other){
    H=other;
}

void MKalmanFilter::setR(const DynamicMatrix<double>& other){
    R=other;
}

void MKalmanFilter::setF(const DynamicMatrix<double>& other){
    F=other;
}

void MKalmanFilter::setB(const DynamicMatrix<double>& other){
    B=other;
}



std::pair<DynamicMatrix<double>, DynamicMatrix<double>> MKalmanFilter::predict(DynamicMatrix<double> _u, DynamicMatrix<double> _B, DynamicMatrix<double> _F, DynamicMatrix<double> _Q){

    // _u, _B, _F, _Q
    if( _B.cols()==1 && _B.rows()==1 && _B(0,0)==0 ){
        _B=B;
    }

    if( _F.cols()==1 && _F.rows()==1 && _F(0,0)==0 ){
        _F=F;
    }
    if( _Q.cols()==1 && _Q.rows()==1 && _Q(0,0)==0){
        _Q=Q;
    }


    // _u, _B, _F, _Q
    if(_B.cols()==1 && _B.rows()==1 && _B(0,0)==0){
      x = _F^x;
    }else{
      x = (_F^x)+(_B^_u);
    }

    // _u, _B, _F, _Q
    if ( _Q.rows()==1 && _Q.cols()==1 && _Q(0,0)==0 ){
        P = ((_F^P)^_F.getTransposed());
    }else{
        P = ((_F^P)^_F.getTransposed()) + _Q;
    }

    //x contains the prior
    return std::make_pair(x,P);
          
}



std::pair<DynamicMatrix<double>, DynamicMatrix<double>> MKalmanFilter::update(DynamicMatrix<double> _z,DynamicMatrix<double> _R, DynamicMatrix<double> _H){
    // error (residual) between measurement and prediction 
    if (_H.cols()==1 && _H.rows()==1 && _H(0,0)==0){
        _H=H;
    }  
    y = _z-(_H^x);  

    //common subexpression S = HPH' + R
    DynamicMatrix<double> PHT = P^_H.getTransposed();  

    // project system uncertainty into measurement space
    if (_R.cols()==1 && _R.rows()==1 && _R(0,0)==0){
        _R=R;
    }
    
    S = (_H^PHT) + _R;

   // map system uncertainty into kalman gain
    K  = PHT^( S.getInverse() );   

    // x = x + Ky
    // predict new x with residual scaled by the kalman gain
    x = x + (K^y);


    DynamicMatrix<double> KH= K^_H;
    assert(KH.rows() == KH.cols());
    DynamicMatrix<double> I_KH = _I-KH;

    auto lhs = (I_KH^P)^I_KH.getTransposed();
    auto rhs =  (K^_R)^K.getTransposed();

    // std::cout << " lhs " << lhs << std::endl;
    // std::cout << " rhs " << rhs << std::endl;
    // std::cout << " P " << P << std::endl;

    // auto __P = lhs + rhs;

    P =  lhs + rhs ;

    return std::make_pair(x,P);   


}