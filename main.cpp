#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <functional>
#include <tuple>
#include "Chapter1/Kalmanfilter.hpp"
#include "Core/matrix.hpp"
#include "Core/constmatrix.hpp"
#include "Core/data.hpp"
#include "Core/dynamicMatrix.hpp"
#include "Chapter2/discreteBayesFilter.hpp"
#include "Chapter2/train.hpp"
#include "Chapter3/probabilities.hpp"
#include "Core/matplotlibcpp.h"
#include "Chapter4/onedimension.hpp"
#include "Chapter4/dogsimulation.hpp"
#include "Chapter5/multivariateGaussians.hpp"


using namespace robotics;
//namespace plt = matplotlibcpp;

int main()
{
    /*******MATRIX OPERATIONS ************************/
    
    Matrix<double, 1, 1> jolin= {{1.0}};
    Matrix<int,1,1> singleMat({{1}});
    Matrix<int, 3, 3> mat2({{2, -3, 1},{2, 0, -1},{1, 4, 5}});
    //mat2+3.0;
    //std::cout << mat2<< std::endl;
    // int det1 = mat2.getDeterminant();
    // auto adjoint  = mat2.getAdjoint();

    // std::cout << " minor " << mat2.getMinor(1,1) << std::endl;
    // std::cout << "determinantYYYY " <<  mat2.getDeterminant() << std::endl;
    // std::cout << " adjugate " << mat2.getAdjugate() << std::endl;
    std::cout << " inverse mat2: " << mat2.getInverse() << std::endl;
    // std::cout << " conver to double "<<  robotics::convert<double>(mat2)<< std::endl;


    /*basketball HEIGHT*/

    Matrix<double,1,2> meanHeight = meanColumnwise(basketballHeight);
    // std::cout<< meanHeight <<std::endl;

    // Matrix<double,1,2> standarDeviation = standardDeviationColumnwise(basketballHeight);
    // std::cout<< standarDeviation<< std::endl;

    /*Usually, measurement errors are distributed normally.
    The Kalman Filter design assumes a normal distribution of the measurement errors.
    Yes, that's correct. Skewness and kurtosis can be used to test for non-normality of data, 
    which is an important assumption for the Kalman filter. 
    If the measurement errors or process noise in a Kalman filter are non-normal, 
    it can affect the accuracy and optimality of the filter's estimates.
    Therefore, before applying the Kalman filter to a dataset, 
    it is important to test for the normality of the measurement errors or process noise. 
    Skewness and kurtosis can be used as indicators of non-normality, 
    along with other statistical tests such as the Shapiro-Wilk test or the Anderson-Darling test.
    If the measurement errors or process noise are significantly non-normal, 
    it may be necessary to use alternative filtering methods or to transform the data to achieve normality. 
    Alternatively, if the non-normality is not severe, 
    it may be possible to modify the Kalman filter to account for the non-normality, 
    such as by using a different likelihood function or by incorporating robust estimators.*/


    //predict using gain guess PAGE 20 book
    double time_step = 1;
    double weight_scale = 0.4; // random number
    double gain_scale = 1.0/3.0; // random number

    double gain_rate = 1; // initial guess
    double weight = 160; // initial guess

    Matrix<double,12,1> estimates;

    for (size_t i = 0; i<weights.getRows();++i){
        //measurement
        const double z = weights.at(i,0);

        // predict new position
        const double prediction = weight + gain_rate*time_step;

        //update filter STATE UPDATE EQUATION kalman filter
        const double innovation = z-prediction; // innovation or residual 
        gain_rate = gain_rate + gain_scale*(innovation/time_step);
        weight = prediction + weight_scale*innovation;

        //append estimates
        estimates.at(i,0)=weight;

    }

    auto estimates2 = g_h_filter(weights,160.0,1.0,0.4,1.0/3.0,1.0);

    std::cout<<" weights "<<weights<<std::endl;
    std::cout<<" estimates "<< estimates<<std::endl;
    std::cout<<" estimates2" << estimates2<<std::endl;

    // GENERATE DATA
    const double initial_value=0.0;
    const double dx_ = 1.0;
    const double noisy_factor = 1.0;
    constexpr size_t count = 30;
    Matrix<double,count,1> Measurements = gen_data<double,count> (initial_value,dx_,noisy_factor);
    std::cout << " measurement "<< Measurements <<std::endl;

    auto estimates3 = g_h_filter(Measurements,initial_value,dx_,0.2,0.02,1.0);
    std::cout << " estimates "<< estimates3 << std::endl;

    // radar tracking
    const double initial_position = 30000; // m 
    const double initial_velocity = 40 ; // m/s 
    //The α−β filter parameters are:
    const double alpha=0.2;
    const double beta=0.1;
    //The track-to-track interval is 5 seconds.
    const double interval = 5;

    auto radarEstimdates = g_h_filter(radarMeasuredPositions,initial_position,initial_velocity,alpha,beta,interval);

    std::cout << " radar estimates "<< radarEstimdates << std::endl;


    // CHAPTER 2 -----------   DISCRETE BAYES FILTER------------------------------------
    std::cout << " pos "<< pos << std::endl;
    std::cout << " halway_map " << hallway_map<<std::endl;

    std::cout << " halfway " << hallway_map.sum()<< std::endl;

    // if sensor gives 1 then 
    std::cout << " dog probabilities "<< std::endl;
    std::cout << calculatePositionProbalities(hallway_map,0.0) <<std::endl;
    std::cout << calculatePositionProbalities(hallway_map,1.0) <<std::endl;



    Matrix<double,1,10> bayesian(0.1);
    auto itsLikelihood = likelihood(hallway_map,1.0,0.75);
    std::cout << update(itsLikelihood,bayesian) << std::endl;


    Matrix<double,1,10> newBay{{.35, .1, .2, .3, 0, 0, 0, 0, 0, .05}};
    std::cout<< perfect_predict(newBay,1) << std::endl;

    Matrix<double,1,10> newBay2{{0., 0., 0., 1., 0., 0., 0., 0., 0., 0.}};
    std::cout<< predict_move(newBay2,2,0.1,0.8,0.1)<<std::endl;
    
    Matrix<double,1,10> newBayesian3{{0, 0, .4, .6, 0, 0, 0, 0, 0, 0}};
    std::cout<< "predict move bayesian3: " << predict_move(newBayesian3,2,0.1,0.8,0.1)<<std::endl;

    Matrix<double,1,3> kernel{{0.1, 0.8, 0.1}};
    std::cout<< "predict move convolution: " << predict_move_convolution(newBayesian3,2.0,kernel) << std::endl;

    Matrix<double,1,10> newBayesian4{{.05, .05, .05, .05, .55, .05, .05, .05, .05, .05}};
    Matrix<double,1,5> kernel2{{.05, .05, .6, .2, .1}};
    std::cout<< predict_move_convolution(newBayesian4,3.0,kernel2)<< std::endl;

    // Measurement Updates -> Predict-> Mesurement Updates -> Predict
    /*
    Initialization
    1. Initialize our belief in the state
    Predict
    1. Based on the system behavior, predict state for the next time step
    2. Adjust belief to account for the uncertainty in prediction
    Update
    1. Get a measurement and associated belief about its accuracy
    2. Compute how likely it is the measurement matches each state
    3. Update state belief with this likelihood
    */
   
    Matrix<double,1,10> prior(0.1);
    constexpr double sensorProb = 0.75;
    double sensorPosition = 1.0; // MEASUREMENT

    // UPDATE SECTION
    auto itsLikelihood_1 = likelihood(hallway_map,sensorPosition,sensorProb);
    auto posterior = update(itsLikelihood_1,prior);

    std::cout << " prior " << prior << std::endl;
    std::cout << " posterior "<< posterior << std::endl;
    
    // PREDICT 
    // perform convolution -> calculate prior 
    Matrix<double,1,3> newKernel{{0.1,0.8,0.1}}; 
    constexpr double sensorMotion = 1.0;
    prior = predict_move_convolution(posterior,sensorMotion,newKernel);

    std::cout << " 1 iteratiorion" << std::endl;
    std::cout << " posterior " << posterior << std::endl;
    std::cout << " calculated prior "<< prior << std::endl;

    // perform update
    sensorPosition = 1.0;
    itsLikelihood_1 = likelihood(hallway_map,sensorPosition,sensorProb);
    posterior = update(itsLikelihood_1,prior);

    /*TRAIN PROBLEM */
    constexpr size_t kernelColumns = 1;
    double sensorAccuracy = .9999;
    Matrix<double,1,kernelColumns> train_kernel = {{1.0}};
    constexpr uint8_t train_move_distance = 4;
    constexpr size_t train_iterations =4;
    Train<kernelColumns> itstrain(10.0,train_kernel,sensorAccuracy);    
    //train_filter(train_iterations,train_kernel,sensorAccuracy,train_move_distance,true);

    // use new kernel
    sensorAccuracy=0.9;
    train_filter(train_iterations,newKernel,sensorAccuracy,train_move_distance,true);

 

    // CHAPTER 3 ------------------------------------------------------
    constexpr size_t numberOfElements = 1000;
    auto mGR = generateRandomNumbers<numberOfElements>();
    //std::cout << mGR << std::endl;

    // EXPECTED VALUE   p_i*X_i for i 1 to N
    double total = 0;
    for(size_t idx=0; idx<numberOfElements;++idx){
        if(mGR.at(0,idx)<=0.8){ // 0.8 probability of being 1
            total+=1;
        }else if(mGR.at(0,idx)<=0.95){ // .15 probability of being 3
            total+=3;
        }else{ // 0.05 probability of being 5
            total+=5;
        }
    }
    std::cout << total/numberOfElements<<std::endl;
    
    // mean and variance
    constexpr ConstMatrix<double, 3, 3> m{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    constexpr double sum1 = MatrixSumUnroll(m);
    std::cout << " sum matrix  " << sum1 << std::endl;
    constexpr double itsmean = mean(m);
    std::cout << " mean "<< itsmean << std::endl;
    constexpr ConstMatrix<double,1,5> X{{1.8, 2.0, 1.7, 1.9, 1.6}};
    constexpr ConstMatrix<double,1,5> Y{{2.2, 1.5, 2.3, 1.7, 1.3}};
    constexpr ConstMatrix<double,1,5> Z{{1.8, 1.8, 1.8, 1.8, 1.8}};
    constexpr double Xmean = mean(X);
    constexpr double Ymean = mean(Y);
    constexpr double Zsmean = mean(Z);
    std::cout << "X mean "<< Xmean << std::endl;
    constexpr double itsSquare = square(X,Xmean);
    constexpr double Xvariance = variance(X);
    std::cout << " X variance "<< Xvariance << std::endl;
    constexpr double XstdDev = standardDeviation(X);
    std::cout << " x stddev "<< XstdDev << std::endl;

    // gaussian 
     const double lower_limit =21.5;
    const double upper_limit =22.5;
    const double mu=22;
    const double sigma=2;

    const double cdf = gaussianIntegral(lower_limit,upper_limit,mu,sigma);
    std::cout << " cdf " << cdf << std::endl;

    // putting all together
    Matrix<double,1,10> itsPrior{{4, 2, 0, 7, 2, 12, 35, 20, 3, 2}};
    itsPrior = itsPrior/itsPrior.sum();
    Matrix<double,1,10> likelihood{{3, 4, 1, 4, 2, 38, 20, 18, 1, 16}};
    likelihood=likelihood/likelihood.sum();

    auto itsPosterior = update(likelihood,itsPrior);
    std::cout << "update " << update(likelihood,itsPrior) << std::endl;
    std::cout << "itsPosterior " << itsPosterior<< std::endl;

    auto meanVar = mean_var(itsPosterior);

    std::cout << "mean " << meanVar.first << " variance " << meanVar.second << std::endl;

    // date age preference 
    // votes for data
    Matrix<double,1,10> dateAge{{0,0,13,43,38,14,8,0,0,0}};
    pollAnalysisAge(dateAge);
    std::cout << " tiktok ad" << std::endl;
    Matrix<double,1,10> tiktokAd{{0,12,32,21,15,15,5,0,0,0}};  
    pollAnalysisAge(tiktokAd); 


    // CHAPTER 4   -----------------------------------------

    double itsMean = 10.0;
    double variance = 1.0;
    std::vector<double> xlim{4,16};
    //plot_gaussian_pdf(itsMean,variance,xlim);
    //plt::show();

    constexpr size_t nElements = 500;
    //auto xs = range(nElements);
    auto ys = randn(nElements)*0.1+10;
    // plt::plot(xs,ys);
    // plt::show();
    std::cout << mean(ys)<< std::endl;

    auto g1 = gaussianS(10, pow(0.2,2));
    auto g2 = gaussianS(15, pow(0.7,2));
    auto g3 = predict(g1,g2);

    std::cout << g3 << std::endl;

    std::cout << g1*g2<<std::endl;

    auto predicted_pos = gaussianS(10., pow(.2,2));
    auto measured_pos = gaussianS(11., pow(.1,2));
    auto estimated_pos = updateODKF(predicted_pos, measured_pos);

    std::cout << estimated_pos << std::endl;

    // 
    auto gn = gaussianS(10, 1);
    auto gnxn = gn*gn;

    // plot_gaussian_pdf(gn.mean,gn.var,xlim);
    // plot_gaussian_pdf(gnxn.mean,gnxn.var,xlim);
    // plt::show();

//N (10.2, 1)×N (9.7, 1).
    auto gLeft=gaussianS(10.2,1);
    auto gRight=gaussianS(9.7,1);
    //plot_products(gLeft,gRight,xlim);

    // FIRST KALMAN FILTER ----------------------------------------
    auto x = gaussianS(0,pow(20,2)); // dogposition HIGH VARIANCE
    double process_var = 2; // 1 
    double velocity = 1;
    double dt = 1; // timestep in seconds
    auto processModel= gaussianS(velocity*dt,process_var); // displacement to add to x
    double sensor_var = pow(300,2); // 2 

    // simulate dog and get measurements
    //auto dog = Dogsimulation(x.mean,velocity,)

    Dogsimulation dog(x.mean,velocity,sensor_var,process_var);
    
    //create list of measurements
    std::vector<double> zs,xs,predictions;
    constexpr size_t nELE = 25;
    for(size_t i=0;i<nELE;++i){
        zs.emplace_back(dog.move_and_sense());
    }

    // plt::bar(zs);
    // plt::plot(zs);
    // plt::show();
    for(const auto& z:zs){
        // predict
        const auto prior = predict(x,processModel);
        // update
        const auto likelihood = gaussianS(z, sensor_var);
        x = updateODKF(likelihood,prior);

        xs.push_back(x.mean);
        predictions.push_back(prior.mean);

        std::cout << " predict "<< prior.mean << " \tvar "<< prior.var ;
        std::cout << "\tmeasurement "<< z <<  "\tupdate "<< x.mean << " \tvar "<< x.var << std::endl;

    }

    std::cout << " actual final positon "<< dog.getX()<< std::endl;
    std::cout << " predicted final position"<< x.mean << std::endl;

    auto xaxis = range(nELE);

/*     plt::plot(xaxis,zs,"k-"); // measurement
    plt::plot(xaxis,xs,"r-"); // estimates 
    //plt::plot(xaxis,predictions,"b-"); // predictiosn
    plt::show(); */

    // INTRODUCTION TO DESIGNING A FILTER ------
    double temp_change = 0;
    double voltage_std = .13;
    double voltProcessVar= pow(0.05,2);
    double actual_voltage = 16.3;

    auto voltX = gaussianS(25,1000); // initial state 
    auto voltProcessModel = gaussianS(0,voltProcessVar); // temperature does not change*
    constexpr size_t voltELE=50;
    
    std::vector<double> voltZS,voltXS,voltPredictions, voltPS;
    voltXS.reserve(voltELE);
    voltPredictions.reserve(voltELE);
    voltPS.reserve(voltELE);

    for(size_t i=0;i<voltELE;++i){
        voltZS.emplace_back(volt(actual_voltage,voltage_std));
    }

    for(const auto& z:voltZS){
        const auto prior = predict(voltX,voltProcessModel);
                // update
        const auto likelihood = gaussianS(z, pow(voltage_std,2));
        voltX = updateODKF(likelihood,prior);

        voltXS.push_back(voltX.mean);
        voltPredictions.push_back(prior.mean);
        voltPS.push_back(voltX.var);

        std::cout << " predict "<< prior.mean << " \tvar "<< prior.var ;
        std::cout << "\tmeasurement "<< z <<  "\tupdate "<< voltX.mean << " \tvar "<< voltX.var << std::endl;


    }

    // plt::plot(voltZS,"k-"); // measurement
    // plt::plot(voltXS,"r-"); // estimates 
    // //plt::plot(voltPredictions,"b-"); // predictiosn
    // plt::show();

    // plt::plot(voltPS);
    // plt::show();

    std::vector<double> W{70.1, 91.2, 59.5, 93.2, 53.5};
    std::vector<double> H{1.8, 2.0, 1.7, 1.9, 1.6};

    std::cout<< covariance(H,W)<<std::endl;


    DynamicMatrix cov({{8.0,0.0},{0.0,3.0}});

    std::cout << " determinant cov " << cov.getDeterminant()<< std::endl;
    std::cout << " inverse cov "<< cov.getInverse() << std::endl;


    std::vector<double> Mu{2., 7.};
    std::vector<double> X1{2.5, 7.3};

    std::cout << " x1-mu "<< X1-Mu << std::endl;

    std::cout << mutlivariateGaussian(X1,Mu,cov) << std::endl;


    double a = 1.0;
    double b = -12.0;
    double c = 31.0;
    auto [root1,root2] = solve_quadratic_equation(a, b, c);

    //std::cout << " root1: "<< root1 << " root2: " << root2 << std::endl;


    //DynamicMatrix<double> covariance_matrix = {{10.0, -2.5}, {-2.5, 8.0}};
    DynamicMatrix<double> covariance_matrix = {{2.0, 0}, {0, 20.0}};
    //DynamicMatrix<double> covariance_matrix = {{2.0, 1.2}, {1.2, 2.0}};
    auto [eigenvectors,eigenvalues] = calculate_eigenvectors(covariance_matrix);
    std::cout  << eigenvectors << std::endl;
    std::cout << eigenvalues << std::endl;
    double v1_x = eigenvectors(0,0);
    double v1_y = eigenvectors(0,1);
    double v2_x = eigenvectors(1,0);
    double v2_y = eigenvectors(1,1);
    plot_ellipse(eigenvalues.at(0), eigenvalues.at(1), v1_x, v1_y, v2_x, v2_y, 1);
 

    return 0;


}