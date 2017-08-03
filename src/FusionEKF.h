#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
    * Constructor.
    */
  FusionEKF();

  /**
    * Destructor.
    */
  virtual ~FusionEKF();

  /**
    * Run the whole flow of the Kalman Filter from here.
    */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
    * Kalman Filter update and prediction math lives in here.
    */
  KalmanFilter ekf_;

  Tools tools;
private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  Eigen::VectorXd x_init_;
  Eigen::MatrixXd F_init_;
  Eigen::MatrixXd P_init_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd H_radar_;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd Q_init_;
};

#endif /* FusionEKF_H_ */
