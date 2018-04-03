#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  bool is_initialized_;

  bool use_laser_;

  bool use_radar_;

  VectorXd x_;

  MatrixXd P_;

  MatrixXd Xsig_pred_;

  long long time_us_;
  double previous_timestamp_;

  double std_a_;

  double std_yawdd_;

  double std_laspx_;

  double std_laspy_;

  double std_radr_;

  double std_radphi_;

  double std_radrd_ ;

  VectorXd weights_;

  int n_x_;

  int n_aug_;

  double lambda_;

  VectorXd x_aug;
    
  MatrixXd P_aug;
    
  MatrixXd Xsig_aug;

  double NIS;
 
  UKF();
 
  virtual ~UKF();

  void ProcessMeasurement(MeasurementPackage meas_package);

  void Prediction(double delta_t);

  void UpdateLidar(MeasurementPackage meas_package);
 
  void UpdateRadar(MeasurementPackage meas_package);
};
#endif