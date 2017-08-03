#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;
  /**
    * Initializes Kalman filter Variables
    * x_ Initial state
    * P_ Initial state covariance
    * F_ Transition matrix
    * H_ Measurement matrix
    * R_ Measurement covariance matrix
    * Q_ Process covariance matrix
    */

  // initializing matrices
  F_init_ = MatrixXd(4, 4);
  P_init_ = MatrixXd(4, 4);
  H_laser_ = MatrixXd(2, 4);
  H_radar_ = MatrixXd(3, 4);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  Q_init_ = MatrixXd(4, 4);
  
  //Initial state transitiion function matrix
  F_init_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //Initial state covariance
  P_init_ << 1, 0, 0,     0,
             0, 1, 0,     0,
             0, 0, 1000,  0,
             0, 0, 0,     1000; 

  //measurement matrix - laser
  //only px and py are accuired by laser measurements
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //measurement matrix - radar
  //r is based on px,py
  //w is based on px,py
  //r' is based on px,py,vx,vy
  H_radar_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  //process covariance matrix
  Q_init_ << 0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0;

}

/**
  * Destructor.
  */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
      */
    // first measurement
    cout << "EKF: " << endl;
    x_init_ = VectorXd(4);
    x_init_ << 0.5, 0.5, 0.1, 0.1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
        * Convert radar from polar to cartesian coordinates and initialize state.
        */
      float r = measurement_pack.raw_measurements_(0);
      float w = measurement_pack.raw_measurements_(1);
      float rdot = measurement_pack.raw_measurements_(2);
      x_init_(0) = r * cos(w);
      x_init_(1) = r * sin(w);
      x_init_(2) = 0;
      x_init_(3) = 0;
      ekf_.Init(x_init_,P_init_,F_init_,H_radar_,R_radar_,Q_init_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
        *Initialize state.
        */
      x_init_(0) = measurement_pack.raw_measurements_(0);
      x_init_(1) = measurement_pack.raw_measurements_(1);
      ekf_.Init(x_init_,P_init_,F_init_,H_laser_,R_laser_,Q_init_);
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  /**
    * Update the state transition matrix F according to the new elapsed time.
    - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

  //set the acceleration noise components
	float noise_ax = 9;
	float noise_ay = 9;
  ekf_.PredictFQCalculation(dt,noise_ax,noise_ay);
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
    * Use the sensor type to perform the update step.
    * Update the state and covariance matrices.
    */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // Calculate the Jacobian of for h(x)
    H_radar_ = tools.CalculateJacobian(ekf_.GetX());
    ekf_.ChangeHR(H_radar_,R_radar_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.ChangeHR(H_laser_,R_laser_);
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
  cout << "x_ = " << ekf_.GetX() << endl;
  cout << "P_ = " << ekf_.GetP() << endl;
}
