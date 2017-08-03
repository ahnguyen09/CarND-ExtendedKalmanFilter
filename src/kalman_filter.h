#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanFilter { 
private:
  // previous state vector
  VectorXd x_prev;

  // state vector
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // state transition matrix
  MatrixXd F_;

  // process covariance matrix
  MatrixXd Q_;

  // measurement matrix
  MatrixXd H_;

  // measurement covariance matrix
  MatrixXd R_;

  // Identity matrix for x_
  MatrixXd I;
public:
  /**
    * Constructor
    */
  KalmanFilter();

  /**
    * Destructor
    */
  virtual ~KalmanFilter();

  /**
    * Init Initializes Kalman filter
    * @param x_in Initial state
    * @param P_in Initial state covariance
    * @param F_in Transition matrix
    * @param H_in Measurement matrix
    * @param R_in Measurement covariance matrix
    * @param Q_in Process covariance matrix
    */
  void Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
      MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in);

  /**
    * Prediction Predicts the state and the state covariance
    * using the process model
    * @param delta_T Time between k and k+1 in s
    */
  void Predict();

  /**
    * Updates the state by using standard Kalman Filter equations
    * @param z The measurement at k+1
    */
  void Update(const VectorXd &z);

  /**
    * Updates the state by using Extended Kalman Filter equations
    * @param z The measurement at k+1
    */
  void UpdateEKF(const VectorXd &z);

  /**
    * Update the state transition matrix F according to the new elapsed time.
    - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
  void PredictFQCalculation(float dt, float noise_ax, float noise_ay);

  /**
    * Convert Cartesian to Polar coordinate for Radar measurment
    */
  VectorXd CartesianToPolar(const VectorXd &x_state); 

  /**
    * Change H and R matrix based on either laser or radar data
    */
  void ChangeHR(MatrixXd &H_in, MatrixXd &R_in);

  /**
    * Public accessor to get state x. Read only
    */
  VectorXd GetX();

  /**
    * Public accessor to get state covariance matrix. Read only
    */
  MatrixXd GetP();

};

#endif /* KALMAN_FILTER_H_ */
