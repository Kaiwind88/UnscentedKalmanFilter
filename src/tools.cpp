#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 
    VectorXd RMSE(4);
    RMSE<< 0,0,0,0;
 
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return RMSE;
    }
    
    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        RMSE += residual;
    }
    
    //calculate the mean
    RMSE = RMSE/estimations.size();
    
    //calculate the squared root
    RMSE = RMSE.array().sqrt();
    
    //return the result
    return RMSE;

}
