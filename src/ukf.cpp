#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    is_initialized_ = false;
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 1;//1.5;//30;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.3;// 0.8;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    n_x_ = 5;
    n_aug_ = 7;
    
    MatrixXd P_ = MatrixXd::Zero(5, 5);
    
    //create augmented mean vector
    x_aug = VectorXd(7);

    //create augmented state covariance
    P_aug = MatrixXd::Zero(7, 7);
    
    //create sigma point matrix
    Xsig_aug = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
    x_ << 0, 0, 0, 0, 0;

    previous_timestamp_ = 0;
    
    Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

    weights_ = VectorXd::Zero(2 * n_aug_ + 1);
    NIS = 0.0;
    
}

UKF::~UKF() {}

/**
 * The latest measurement data of either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    /*  Initialization */
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
        (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {
        if (!is_initialized_) {
        /**
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         */
        // first measurement
        // 0 for lidar, 1 for radar
        P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 10, 0, 0,
        0, 0, 0, 50, 0,
        0, 0, 0, 0, 3;

        if (meas_package.sensor_type_ == 0) {
            x_(0) =meas_package.raw_measurements_(0);
            x_(1) =meas_package.raw_measurements_(1);
        } else if (meas_package.sensor_type_ == 1) {
            float px_radar = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
            float py_radar = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
            float v_abs_radar = meas_package.raw_measurements_(2);
            x_(0) = px_radar;
            x_(1) = py_radar;
        }
        
        previous_timestamp_ =  meas_package.timestamp_ ;
        cout<< "Initialized done!";
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
        double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt in seconds
        previous_timestamp_ =  meas_package.timestamp_ ;
        Prediction(dt);
        
        if (meas_package.sensor_type_ == 0) {
            UpdateLidar(meas_package);
        } else if (meas_package.sensor_type_ == 1) {
            UpdateRadar(meas_package);
        

    }
    // RADAR: 1
    // LASER: 0
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {
  /**
    Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
  **/
  /* UKF Augumentated Generate Sigma Point */  
    //create augmented mean state
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    lambda_ = 3 - n_aug_;
    
    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    //cout << "P_aug"<<P_aug<<endl;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    //cout << "sig point estimate"<<Xsig_aug <<endl;
    
    /* Sigma Point Prediction */
    //create matrix with predicted sigma points as columns
    
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
    //cout << "sig point predict"<<Xsig_pred_<<endl;
    
    /* Calculate Predicted Mean x_ and Covariance P_*/
    
    // set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
        double weight = 0.5/(n_aug_+lambda_);
        weights_(i) = weight;
    }
    //cout << "weights"<<weights_<<endl;
    //predicted state mean
    x_.fill(0.0);
    //cout << "x_ is 0 now: "<< x_ <<endl;

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    }
    
    //x_ << x_aug.head(5);
    //cout << "x_ after predict step:"<< x_ <<endl;

    //predicted state covariance matrix
    P_.fill(0.0);
    //P_aug.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //cout <<"Get short column" <<sig_7.head(5)<<endl;

        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
        //cout <<"P_" <<P_<<endl;
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Modify the state vector x_, and covariance P_

  Use lidar data to update the belief about the object's position. Modify the state vector, x_, and covariance, P_.

  **/
    
    //set measurement dimension, lidar can measure px, py
    int n_z = 2;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
       
        // measurement model
        Zsig(0,i) = p_x;                        
        Zsig(1,i) = p_y;                                 
        
    }
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;
    S = S + R;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
   
    // * Radar update
     
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
                // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    NIS = z_diff.transpose() * S.inverse()*z_diff;
    cout << "LIDAR NIS: " << NIS<<endl;
    
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    
    
     /* Modify the state vector x_, and covariance P_ */
     
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;
    S = S + R;
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    
    /* Radar update */
 
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    NIS = z_diff.transpose() * S.inverse()*z_diff;
    cout << "RADAR NIS: " << NIS<<endl;
}
