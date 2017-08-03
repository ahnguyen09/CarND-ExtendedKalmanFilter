#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

/**
  * Takes in ground truth and estimation state (px,py,vx,vy) and calculates root mean squred error
  */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
                    
  VectorXd rmse(4);
	VectorXd err(4);
	rmse << 0,0,0,0;

  /**
  * check the validity of the following inputs:
  * the estimation vector size should not be zero
  * the estimation vector size should equal ground truth vector size
  */
	if (estimations.size() != 0 and estimations.size() == ground_truth.size()) {
    	for(int i=0; i < estimations.size(); ++i){
            err = estimations[i] - ground_truth[i];
            err = err.array()*err.array();
            //accumulate squared residuals
            rmse = rmse + err; 
    	}
    
    	//calculate the mean
    	rmse = rmse/estimations.size();
    
    	//calculate the squared root
    	rmse = rmse.array().sqrt();	
	}
	else {
	    cout << "Mismatch size or empty" << endl;
	}

	//return the result
	return rmse;
}

/**
  * Used in measurment update step for S,K,P calculations
  * Takes in Cartesian x state (px,py,vx,vy) and return polar Jacobian matrix 3x4
  */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  //Size of Jacobian is 3x4 because we are mapping px,py,vx,vy to r,w,r'; 3x4 * 4x1 -> 3x1 
  MatrixXd Hjb(3,4);
  //recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  
  float px2_py2; //defined as px^2 + py^2
  float sqrt_px2_py2; //defined as sqrt(px^2 + py^2)
  float px2_py2_32; //defined as (px^2 + py^2)^(3/2)

  //check divide by 0
  if (px == 0 and py == 0) {
    cout << "CalculateJacobian () - Error - Divide by 0";
    Hjb << 0,0,0,0,
           0,0,0,0,
           0,0,0,0;
    return Hjb;
  }

  px2_py2 = pow(px,2) + pow(py,2);
  sqrt_px2_py2 = sqrt(px2_py2);
  px2_py2_32 = px2_py2 * sqrt_px2_py2;

  Hjb << px/sqrt_px2_py2,               py/sqrt_px2_py2,                0,                0,
	       -py/px2_py2,                   px/px2_py2,                     0,                0,
	       (py*(vx*py-vy*px))/px2_py2_32, (px*(vy*px-vx*py))/px2_py2_32,  px/sqrt_px2_py2,  py/sqrt_px2_py2;

  return Hjb;
}
