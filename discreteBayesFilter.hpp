#ifndef DISCRETEBAYESFILTER_HPP
#define DISCRETEBAYESFILTER_HPP

#include "matrix.hpp"
#include "Utils.hpp"
#include "train.hpp"

using namespace robotics;

Matrix<double,1,10> pos{{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1}};

Matrix<double,1,10> hallway_map{{1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0}};

/// @brief OTHER is made of ones and zeros. If the sensor is 1 all good. If the sensorValue is zero 
// make ones to zeros and zeros to ones. Normalize before returning
/// @tparam T 
/// @tparam ROWS 
/// @tparam COLS 
/// @param other 
/// @param sensorValue 
/// @return 
template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> calculatePositionProbalities(const Matrix<T,ROWS,COLS>& other, T sensorValue){
    Matrix<T,ROWS,COLS> pos_;
    if(sensorValue==0){
        for(size_t row = 0; row<ROWS;++row){
            for(size_t col=0; col<COLS;++col){
                pos_.at(row, col) = ( (other.at(row, col) == 0) ? 1.0 : 0.0 );
            }
        }        
    }else{
        pos_ = other*1.0;
    }
    return pos_/pos_.sum();
}

/// @brief  POSTERIOR = LIKELIHOOD * PRIOR  / NORMALIZATION
/// @tparam T 
/// @tparam ROWS 
/// @tparam COLS 
/// @param likelihood 
/// @param prior 
/// @return UPDATE VERSION
template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> update(const Matrix<T,ROWS,COLS>& likelihood, const Matrix<T,ROWS,COLS>& prior){
    auto update_belief = likelihood*prior;      
    return update_belief/static_cast<T>(update_belief.sum());
}

/// @brief COMPUTE LIKELIHOOD FOR HALLWAY
    //  Iterate over the map of positions defined by hallway_map and each element that matches 
    //  SensorReading increase it by a scale factor defined by SCALE
/// @tparam T 
/// @tparam ROWS 
/// @tparam COLS 
/// @param hallway_map 
/// @param sensorReading 
/// @param sensorReadingProbability 
/// @return likelihood
template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> likelihood(const Matrix<T,ROWS,COLS>& hallway_map, T sensorReading , T sensorReadingProbability){

    
    T scale = sensorReadingProbability/T(1-sensorReadingProbability);
    Matrix<T,ROWS,COLS> likelihood(1); // at the beginning
        for(size_t row = 0; row<ROWS;++row){
            for(size_t col=0; col<COLS;++col){
                const T& hallway_map_value = hallway_map.at(row,col);
                const T& likelihood_value = likelihood.at(row,col);
                likelihood.at(row,col) = ( hallway_map_value==sensorReading?likelihood_value*scale:likelihood_value );
            }
        }        
    return likelihood;
}

/// @brief move the position by `move` spaces, where positive is to the right, and negative is to the left
/// @tparam T 
/// @tparam ROWS 
/// @tparam COLS 
/// @param belief 
/// @param move 
/// @return 
template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> perfect_predict(const Matrix<T,ROWS,COLS>& belief, size_t move){
    // we are only shifting columwise
    Matrix<T,ROWS,COLS> perfectPredict;
    for(size_t row=0;row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            const size_t indexCol = modulo<int>(int(col)-int(move), int(COLS));
            //std::cout << " col " << col << " move " << move << " indexCol: " << indexCol<< std::endl;
            perfectPredict.at(row,col) = belief.at(row,indexCol);
            //std::cout<<" belief at (" << row << ","<< indexCol<<") "<< belief.at(row,indexCol)<< std::endl;
        }
    }
    return perfectPredict;
}

template<typename T, size_t ROWS, size_t COLS>
Matrix<T,ROWS,COLS> predict_move(const Matrix<T,ROWS,COLS>& belief, size_t move, T prob_under, T prob_correct, T prob_over){
    Matrix<T,ROWS,COLS> predict;
    for(size_t row=0; row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){

            const size_t indexCorrect= modulo(int(col)-int(move),int(COLS));
            const size_t indexUnder=modulo(int(col)-int(move)+1,int(COLS));
            const size_t indexOver=modulo(int(col)-int(move)-1,int(COLS));
            predict.at(row,col) =   belief.at(row,indexCorrect)*prob_correct
                                    +belief.at(row,indexOver)*prob_over
                                    +belief.at(row,indexUnder)*prob_under;
        }
    }

    return predict;
}

// 
// 

/// @brief generalization of predict_move(...) using a KERNEL instead of prob_under,prob_correct,prob_over
    //move -> offset
/// @tparam T 
/// @tparam ROWS 
/// @tparam COLS 
/// @tparam KERNELCOLS 
/// @param belief 
/// @param offset 
/// @param kernel 
/// @return 
template<typename T, size_t ROWS, size_t COLS, size_t KERNELCOLS>
Matrix<T,ROWS,COLS> predict_move_convolution(const Matrix<T,ROWS,COLS>& belief, T offset, const Matrix<T,ROWS,KERNELCOLS>& kernel){
    Matrix<T,ROWS,COLS> predict(0);
    //assert(KERNELCOLS%2==1 && KERNELCOLS>=3 && KERNELCOLS<COLS); // odd number, larger than 3 and smaller than kenerlcols 
    const size_t kernel_halfwidth = (KERNELCOLS-1)/2;

    for(size_t row=0; row<ROWS;++row){
        for(size_t col=0;col<COLS;++col){
            T& elementToBeModified = predict.at(row,col);
            for(size_t kernelcol=0; kernelcol<KERNELCOLS;++kernelcol){
                const int8_t preShift = kernel_halfwidth-kernelcol;
                const size_t adjustedIndex = modulo(int(col)-int(offset)+preShift,int(COLS));

                //index = (i + (width-k) - offset) % N
                //std::cout << "col: " << col << " kernelcol: "<< kernelcol << " adjustedIndex "<< adjustedIndex << std::endl;
                elementToBeModified+=belief.at(row,adjustedIndex)*kernel.at(row,kernelcol);
            }


        }
    }
    return predict;
}

/// @brief the comand_distance is what the robot should have moved 
/// robot.move() simulates a noise of the movement of the robot 
/// the predict_move_convolution does not know this so it calculates the next position with the kernel, an move_distance
/// robot.sense() has also noise incorporated. This is use for the measurement
/// @tparam KERNELCOLS 
/// @param iterations 
/// @param kernel 
/// @param sensor_accuracy 
/// @param move_distance 
/// @param doPrint 
template<size_t KERNELCOLS>
void train_filter(size_t iterations, Matrix<double,1,KERNELCOLS> kernel, double sensor_accuracy, int move_distance, bool doPrint=true){
    // track = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    //          prior = np.array([.9] + [0.01]*9)
    
    constexpr uint8_t track_len = 10;
    Matrix<double,1,track_len> track{{0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9.}};
    Matrix<double,1,track_len> prior{{0.9,.01,.01,.01,.01,.01,.01,.01,.01,.01}};
    auto posterior = prior;
    prior = prior/prior.sum(); // normalize it 

    Train<KERNELCOLS> robot(track_len,kernel,sensor_accuracy); 

    //for i in range(iterations):
    for(size_t iter = 0; iter<iterations;++iter){
        // the command "move_distance "
        robot.move(move_distance); // update pos member variable from robot which takes NOISE into account
        // this value is used in  SENSE since 
        
        // perform predict
        //double itsMoveDistance = move_distance;
        prior = predict_move_convolution(posterior,static_cast<double>(move_distance),kernel);

        // call sense from robot which also takes NOISE into account
        const double measuredPosition = robot.sense();

        const Matrix<double,1,track_len> itsLikeliehood = likelihood(track,measuredPosition,sensor_accuracy);
        posterior = update(itsLikeliehood,prior);
        auto pair = getMaxElementIndex(posterior); // first ROW second COL 

        if(doPrint){

            const size_t robot_pos = robot.getTrainPos();
            std::cout << " iteration " << iter << " at position " << robot_pos << " sensed " 
            << measuredPosition << "  at position "<< track.at(0,robot_pos)<< std::endl;

            std::cout << "          estimated position is " << pair.second << " with confidence "<< 
            posterior.at(pair.first,pair.second)*100 << "% " << std::endl;
        
        }
    }

}

#endif 