#include "kalman_filter.h"
#include <math.h>
#include <iostream>

const float PI_2 = 2*M_PI;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

  I = MatrixXd::Identity(4, 4);
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_; // 2x1 = 2x4 * 4x1
	VectorXd y = z - z_pred; // 2x1
	
	MatrixXd Ht = H_.transpose(); 
  MatrixXd PHt = P_ * Ht; // 4x2 =  4x4 * 4x2
	MatrixXd S = H_ * PHt + R_; // 2x2 = 2x4 * 4x2 + 2x2
	MatrixXd Si = S.inverse(); 
	MatrixXd K = PHt * Si; // 4x2 = 4x2 * 2x2

  x_ = x_ + (K * y); // 4x1 = 4x2 * 2x1 
	P_ = (I - K * H_) * P_; // 4x4 = (4x4 - (4x2 * 2x4)) * 4x4
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = CartesianToPolar(x_); // 3x1	
	VectorXd y = z - z_pred;

	//Normalize angle between -pi and pi
	int tquot = y(1)/M_PI;
	y(1) = y(1) - (tquot * M_PI);
  y(1) = y(1) - (tquot%2)*M_PI;

	MatrixXd Ht = H_.transpose(); 
  MatrixXd PHt = P_ * Ht; // 4x3 = 4x4 * 4x3
	MatrixXd S = H_ * PHt + R_; // 3x3 = 3x4*4x3 + 3x3
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si; // 4x3 = 4x3 * 3x3

  x_ = x_ + (K * y); // 4x1 = 4x3 * 3x1
	P_ = (I - K * H_) * P_; // 4x4 = (4x4 - (4x3 * 3x4)) * 4x4
}

VectorXd KalmanFilter::CartesianToPolar(const VectorXd &x_state){
	VectorXd z(3);
	float px, py, vx, vy;
  px = x_state[0];
  py = x_state[1];
  vx = x_state[2];
  vy = x_state[3];

	float r,w,rdot;
	float px2 = pow(px,2);
	float py2 = pow(py,2);
  float sqrt_px2py2 = sqrt(px2 + py2);

	r = sqrt_px2py2;
	w = atan2(py,px); // result is between -pi and pi 
	if (r < 0.000001) r = 0.000001; // prevent rdot calculation: divide by 0 error
	rdot = (px*vx + py*vy)/r;
	z<< r,w,rdot;
	return z;
}

void KalmanFilter::PredictFQCalculation(float dt, float noise_ax, float noise_ay) {
  float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	F_(0, 2) = dt;
	F_(1, 3) = dt;

	//set the process covariance matrix Q
	Q_ <<  (dt_4*noise_ax)/4, 0, (dt_3*noise_ax)/2, 0,
			   0, (dt_4*noise_ay)/4, 0, (dt_3*noise_ay)/2,
			   (dt_3*noise_ax)/2, 0, dt_2*noise_ax, 0,
			   0, (dt_3*noise_ay)/2, 0, dt_2*noise_ay;
}

void KalmanFilter::ChangeHR(MatrixXd &H_in, MatrixXd &R_in) {
	H_ = H_in;
	R_ = R_in;
}

VectorXd KalmanFilter::GetX() {return x_;}

MatrixXd KalmanFilter::GetP() {return P_;}
