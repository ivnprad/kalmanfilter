#ifndef MKALMANFILTER_HPP
#define MKALMANFILTER_HPP

#include "../Core/dynamicMatrix.hpp"

//class DynamicMatrix;

class MKalmanFilter{
    public:
        MKalmanFilter( int _dim_x, int _dim_z, int _dim_u=0);
        std::pair<DynamicMatrix<double>, DynamicMatrix<double>> predict(
            DynamicMatrix<double> u=DynamicMatrix<double> (1,1), 
            DynamicMatrix<double> B=DynamicMatrix<double> (1,1),
            DynamicMatrix<double> F=DynamicMatrix<double> (1,1),
            DynamicMatrix<double> Q=DynamicMatrix<double> (1,1)
        );
        std::pair<DynamicMatrix<double>, DynamicMatrix<double>> update(
            DynamicMatrix<double> _z,
            DynamicMatrix<double> _R=DynamicMatrix<double> (1,1),
            DynamicMatrix<double> _H=DynamicMatrix<double> (1,1)
        );

        int getDimX() const;
        int getDimZ() const;

        DynamicMatrix<double> getX() const;
        DynamicMatrix<double> getP() const;
        DynamicMatrix<double> getX_PRIOR() const;
        DynamicMatrix<double> getP_PRIOR() const;
        DynamicMatrix<double> getX_POST() const;
        DynamicMatrix<double> getP_POST() const;
        DynamicMatrix<double> getQ() const;
        DynamicMatrix<double> getR() const;
        DynamicMatrix<double> getF() const;
        DynamicMatrix<double> getH() const;
        DynamicMatrix<double> getK() const;
        DynamicMatrix<double> getY() const;
        DynamicMatrix<double> getS() const;
        DynamicMatrix<double> getSI() const;
        DynamicMatrix<double> getM() const;
        DynamicMatrix<double> getB() const;
        DynamicMatrix<double> getZ() const;


        double getLogLikeihood() const;
        double getLikeihood() const;
        double getMahalanobis() const;
        double getAlpha() const;

       
        void setX(const DynamicMatrix<double>& other);
        void setP(const DynamicMatrix<double>& other);
        void setQ(const DynamicMatrix<double>& other);
        void setH(const DynamicMatrix<double>& other);
        void setR(const DynamicMatrix<double>& other);
        void setF(const DynamicMatrix<double>& other);
        void setB(const DynamicMatrix<double>& other);


    private:

        int dim_x;
         /*
        Number of state variables for the Kalman filter. For example, if
        you are tracking the position and velocity of an object in two
        dimensions, dim_x would be 4.
        This is used to set the default size of P, Q, and u
        */       
        int dim_z;
        /*
        Number of of measurement inputs. For example, if the sensor
        provides you with position in (x,y), dim_z would be 2.        
        */
        int dim_u;
        /*
        size of the control input, if it is being used.
        Default value of 0 indicates it is not used.        
        */

        DynamicMatrix<double> x; //Current state estimate. Any call to update() or predict() updates this variable.
        DynamicMatrix<double> P; // uncertainty covariance Current state covariance matrix. Any call to update() or predict() updates this variable.   
        DynamicMatrix<double> Q; // process covariance 
        DynamicMatrix<double> B; // control transition matrix 
        DynamicMatrix<double> F; // state transition matrix 
        DynamicMatrix<double> H; // measurement function
        DynamicMatrix<double> R; // state uncertainty
        double _alpha_sq;       // fading memory control 
        DynamicMatrix<double> M; // process-measurement cross correlation
        DynamicMatrix<double> z; 

        /*
        gain and residual are computed during the innovation step. We
        save them so that in case you want to inspect them for various
        purposes*/
        DynamicMatrix<double> K; // kalman Gain
        DynamicMatrix<double> y; // 
        DynamicMatrix<double> S; // system uncertainty
        DynamicMatrix<double> SI; // inverse sytem uncertainty

        const DynamicMatrix<double> _I; // # identity matrix. Do not alter this.

        //# these will always be a copy of x,P after predict() is called
        DynamicMatrix<double> x_prior;
        DynamicMatrix<double> P_prior;

        // these will always be a copy of x,P after update() is called
        DynamicMatrix<double> x_post;
        DynamicMatrix<double> P_post;

        //Only computed only if requested via property
        double log_likelihood;
        double likelihood;
        double _mahalanobis;
       


/*

Predict

Univariate
μ¯=μ+μfx
σ¯2=σ2x+σ2fx

Univariate(Kalman form)
x¯=x+dx
P¯=P+Q

Multivariate

x¯=Fx+Bu
P¯=FPF𝖳+Q

Without worrying about the specifics of the linear algebra, we can see that:

x,P
 are the state mean and covariance. They correspond to x
 and σ2
.

F
 is the state transition function. When multiplied by x
 it computes the prior.

Q
 is the process covariance. It corresponds to σ2fx
.

B
 and u
 are new to us. They let us model control inputs to the system

*/
    

};

#endif 