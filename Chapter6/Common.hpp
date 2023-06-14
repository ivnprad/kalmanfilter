#ifndef COMMON_HPP
#define COMMON_HPP

#include <vector>
#include <cmath>
#include "../Core/Utils.hpp"
#include "../Core/dynamicMatrix.hpp"
#include "MKalmanfilter.hpp"
#include "../Core/matplotlibcpp.h"

template<typename T=double>
std::pair<std::vector<T>,std::vector<T>> compute_dog_data( const T z_var, const T process_var, const int count=1,const T dt=1.){
   T x =0.;
   const T vel=1.0;
   const T z_std = sqrt(z_var);
   const T p_std = sqrt(process_var);
   std::vector<T> xs,zs;
   for(int idx=0; idx<count;++idx){

      const T v = vel +  ( myRandom(-1.0,1.0)*p_std);
      x+= v*dt;
      xs.emplace_back(x);
      zs.emplace_back(x+ myRandom(-1.0,1.0)*z_std); 

   }

   return std::make_pair(xs,zs);

}

// def run(x0=(0.,0.), P=500, R=0, Q=0, dt=1.0, 
//         track=None, zs=None,
//         count=0, do_plot=True, **kwargs):

template<typename T=double>
std::pair<std::vector<DynamicMatrix<T>>,std::vector<DynamicMatrix<T>>> run(DynamicMatrix<T> x, DynamicMatrix<T> P, DynamicMatrix<T> R, DynamicMatrix<T> Q, 
      double dt=1.0, std::vector<T> xs = std::vector<double> (0), std::vector<T> zs = std::vector<double> (0), int count=0){

   std::cout << " what the hell " << std::endl;
   //Simulate dog if no data provided.    
   auto [xs_,zs_] = compute_dog_data(R(0,0),Q(0,0),count);
   xs = xs_;
   zs = zs_;
   std::vector<T> state_position, pos_var, vel_var;
   namespace plt = matplotlibcpp;


   //create the kalman filter
   MKalmanFilter kf = pos_vel_filter(x,P,R,Q,dt);

   std::vector<DynamicMatrix<T>> _xs;
   std::vector<DynamicMatrix<T>> COV;
   
   // run the kalman filter and store the results
   for( const auto& element:zs){
      kf.predict();
      kf.update(DynamicMatrix<double> ({{element}}));
      _xs.push_back(kf.getX());
      state_position.push_back(kf.getX()(0,0));
      COV.push_back(kf.getP());
      pos_var.push_back(kf.getP()(0,0));
      vel_var.push_back(kf.getP()(1,1));
   }

   plt::plot(xs); //track
   plt::plot(zs); //measurement
   plt::plot(state_position); //kalman filtered
   plt::plot(pos_var);
   plt::plot(vel_var);
   plt::show();

   return std::make_pair(_xs,COV);

}



template<typename T=double>
std::pair<DynamicMatrix<T>, DynamicMatrix<T>> predictMKF(DynamicMatrix<T> x, DynamicMatrix<T> P, DynamicMatrix<T> F, DynamicMatrix<T> Q , DynamicMatrix<T> B = DynamicMatrix(1,1), DynamicMatrix<T> u = DynamicMatrix(1,1)){
   
   assert(F.rows()>=2);
   assert(F.cols()>=2);

   if(B.rows()==1 && B.cols()==1 ){
      x = F^x;
   }else{
      x = (F^x)+(B^u);
   }
   
   if ( Q.rows()==1 && Q(0,0)==0 ){
      P = ((F^P)^F.getTransposed());
   }else{
      P = ((F^P)^F.getTransposed()) + Q;
   }

   return std::make_pair(x,P);

}



template<typename T=double>
std::pair<DynamicMatrix<T>, DynamicMatrix<T>> updateMKF(DynamicMatrix<T> x, DynamicMatrix<T> P, DynamicMatrix<T> z, DynamicMatrix<T> R, DynamicMatrix<T> H){
   
   //y = z - Hx
   // error (residual) between measurement and prediction   
   DynamicMatrix<T> y = z-(H^x);
   // std::cout << " y " << y << std::endl;
   // std::cout << " z- H^x "<< (z-(H^x)) << std::endl;

   //common subexpression S = HPH' + R
   DynamicMatrix<T> PHT = P^H.getTransposed();

   // project system uncertainty into measurement space
   DynamicMatrix<T> S = (H^PHT) + R;

   // map system uncertainty into kalman gain
   DynamicMatrix<T> K  = PHT^( S.getInverse() );

   // x = x + Ky
   // predict new x with residual scaled by the kalman gain
   x = x + (K^y);


   DynamicMatrix<T> KH= K^H;
   assert(KH.rows() == KH.cols());
   DynamicMatrix<T> I = DynamicMatrix<double>::IdentityMatrixCreator<>::create(KH.rows());
   DynamicMatrix<T> I_KH = I-KH;

   P =  ( (I_KH^P)^I_KH.getTransposed() )+ ( (K^R)^K.getTransposed() ) ;


   return std::make_pair(x,P); 

}


// z = 1.
// x, P = update(x, P, z, R, H)
// print('x =', x)
// x = [ 1.085 -0.64 ]


template<typename T=double>
DynamicMatrix<T> Q_discrete_white_noise(const size_t dim, const T dt=1.0, const T var=1.0, const size_t block_size=1, bool order_by_dim = true){
    /*
    Returns the Q matrix for the Discrete Constant White Noise
    Model. dim may be either 2, 3, or 4 dt is the time step, and sigma
    is the variance in the noise.

    Q is computed as the G * G^T * variance, where G is the process noise per
    time step. In other words, G = [[.5dt^2][dt]]^T for the constant velocity
    model.

    Parameters
    -----------

    dim : int (2, 3, or 4)
        dimension for Q, where the final dimension is (dim x dim)

    dt : float, default=1.0
        time step in whatever units your filter is using for time. i.e. the
        amount of time between innovations

    var : float, default=1.0
        variance in the noise

    block_size : int >= 1
        If your state variable contains more than one dimension, such as
        a 3d constant velocity model [x x' y y' z z']^T, then Q must be
        a block diagonal matrix.

    order_by_dim : bool, default=True
        Defines ordering of variables in the state vector. `True` orders
        by keeping all derivatives of each dimensions)

        [x x' x'' y y' y'']

        whereas `False` interleaves the dimensions

        [x y z x' y' z' x'' y'' z'']


    Examples
    --------
    >>> # constant velocity model in a 3D world with a 10 Hz update rate
    >>> Q_discrete_white_noise(2, dt=0.1, var=1., block_size=3)
    array([[0.000025, 0.0005  , 0.      , 0.      , 0.      , 0.      ],
           [0.0005  , 0.01    , 0.      , 0.      , 0.      , 0.      ],
           [0.      , 0.      , 0.000025, 0.0005  , 0.      , 0.      ],
           [0.      , 0.      , 0.0005  , 0.01    , 0.      , 0.      ],
           [0.      , 0.      , 0.      , 0.      , 0.000025, 0.0005  ],
           [0.      , 0.      , 0.      , 0.      , 0.0005  , 0.01    ]])

    References
    ----------

    Bar-Shalom. "Estimation with Applications To Tracking and Navigation".
    John Wiley & Sons, 2001. Page 274.
    */

   DynamicMatrix<double> Q(1,1);

   if ( dim <2  || dim>4){
      throw std::runtime_error("dim must be between 2 and 4");
   }

   switch (dim)
   {
      case 2:
         Q = DynamicMatrix({
            {0.25 * pow(dt, 4), 0.5 * pow(dt, 3)}, 
            {0.5 * pow(dt, 3), pow(dt,2)}
            });
         break;

      case 3:
         Q = DynamicMatrix(
            {{0.25 * pow(dt, 4), 0.5 * pow(dt, 3), 0.5 * pow(dt, 2)},
              {0.5 * pow(dt, 3),    pow(dt,2),       dt},
              { .5*pow(dt,2),       dt,        1}}
         );
         break; 
      case 4:
         Q = DynamicMatrix(
            {
               {pow(dt,6)/36, pow(dt,5)/12,  pow(dt,4)/6,   pow(dt,3)/6},
               {pow(dt,5)/12, pow(dt,4)/4,   pow(dt,3)/2,   pow(dt,2)/2},
               {pow(dt,4)/6,  pow(dt,3)/2,   pow(dt,2),     dt},
               {pow(dt,3)/6,  pow(dt,2)/2,   dt,            1.0}
            }
         );        
   
      default:
         break;
   }

   // if (order_by_dim){
   //    Q*=var; // this causes problems?
   // }
   
   return Q*var;


}


MKalmanFilter univariate_filter(const DynamicMatrix<double>& x, const DynamicMatrix<double>& P, const DynamicMatrix<double>& R, const DynamicMatrix<double>& Q= DynamicMatrix(1,1 )){
   const int dim_x=1;
   const int dim_z=1;
   MKalmanFilter kf(dim_x,dim_z);
   kf.setX(x);
   kf.setP(P);
   kf.setQ(Q);
   kf.setR(R);
   kf.setH(DynamicMatrix<double>({{1.0}}));
   kf.setF(DynamicMatrix<double> ({{1.0}}));
   kf.setB(DynamicMatrix<double> ({{1.0}}));

   return kf;
}

MKalmanFilter pos_vel_filter(const DynamicMatrix<double>& x, const DynamicMatrix<double>& P, const DynamicMatrix<double>& R, const DynamicMatrix<double>& Q= DynamicMatrix(1,1), const double dt=1.0 ){
    // Returns a KalmanFilter which implements a
    // constant velocity model for a state [x dx].T
    // 
   const int dim_x=2; // state variables to watch 
   const int dim_z=1; // measurement
   MKalmanFilter kf(dim_x,dim_z);
   kf.setX(x);
   

   //state transition matrix
   DynamicMatrix<double> F({{1.0,dt},{0.,1.0}});
   kf.setF(F);
   DynamicMatrix<double> H({{1.0,0.0}});
   kf.setH(H);

   /*
   KalmanFilter initializes R, P, and Q to the identity matrix, 
   so kf.P *= P is one way to quickly assign all of the diagonal elements to the same scalar value. 
   */
   kf.setR(kf.getR()*R);

   std::cout << " getP " << kf.getP() << std::endl;
   std::cout << " P " << P  << std::endl;

   if(P.cols()==kf.getP().cols() && P.rows()==kf.getP().rows()){
      kf.setP(kf.getP()*P);
   }else if (P.cols()==1&&P.rows()==1)
   {
      kf.setP(kf.getP()*P(0,0));
   }
   
   

   if(Q.rows()==1&&Q.cols()==1){
      auto itsQ = Q_discrete_white_noise(dim_x,dt,Q(0,0));
      kf.setQ(itsQ);

   }else{
      kf.setQ(Q);
   }

   return kf;

}

void print_kf(const MKalmanFilter& KF ){

    std::cout<< " dim_x : " << KF.getDimX() << std::endl;
    std::cout<< " dim_z : " << KF.getDimZ() << std::endl;
    std::cout << " x: " << KF.getX() << std::endl;
    std::cout << " P: " << KF.getP() << std::endl;

    std::cout << " x_prior "<< KF.getX_PRIOR() << std::endl;
    std::cout << " P_prior "<< KF.getP_PRIOR() << std::endl;
    std::cout << " x_post "<< KF.getX_POST()<< std::endl;
    std::cout << " P_post "<< KF.getP_POST()<< std::endl;

    std::cout << " F: " << KF.getF() << std::endl;
    std::cout << " Q: " << KF.getQ() << std::endl;
    std::cout << " R: " << KF.getR() << std::endl;
    std::cout << " H: " << KF.getH() << std::endl;
    std::cout << " K: " << KF.getK() << std::endl;
    std::cout << " y: " << KF.getY() << std::endl;
    std::cout << " S: " << KF.getS() << std::endl;


}


std::ostream& operator<<(std::ostream& stream, const MKalmanFilter& KF){
   // for(size_t i=0;i<ROWS;++i){
   //    for(size_t j=0;j<COLS;++j){
            
   // }

   stream <<"dim_x : " << KF.getDimX() << "\n";
   stream<< "dim_z : " << KF.getDimZ() << "\n";
   stream << "x: " << KF.getX() << "\n";
   stream << "P: " << KF.getP() << "\n";  

   stream<< "x_prior "<< KF.getX_PRIOR() << "\n";
   stream<< "P_prior "<< KF.getP_PRIOR() << "\n";
   stream<< "x_post "<< KF.getX_POST() << "\n";
   stream<< "P_post "<< KF.getP_POST() << "\n";  

   stream << " F: " << KF.getF() << "\n";
   stream << " Q: " << KF.getQ() << "\n";
   stream << " R: " << KF.getR() << "\n";
   stream << " H: " << KF.getH() << "\n";
   stream << " K: " << KF.getK() << "\n";
   stream << " y: " << KF.getY() << "\n";
   stream << " S: " << KF.getS() << "\n";
   stream << " SI: " << KF.getSI() << "\n";
   stream << " M." << KF.getM() << "\n";
   stream << " B." << KF.getB() << "\n";
   stream << " Z." << KF.getZ() << "\n";

   stream << " log-likeihood " << KF.getLogLikeihood() << std::endl;
   stream << " likeihood "<< KF.getLikeihood() << std::endl;
   stream << " mahalanobis "<< KF.getMahalanobis() << std::endl;
   stream << " alpha " << KF.getAlpha() << std::endl;
   

   stream << "\n";

   // }
   return stream;
}



#endif 