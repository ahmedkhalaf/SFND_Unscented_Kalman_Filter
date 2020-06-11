#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

void NormalizeAngle(double& phi)
{
  phi = atan2(sin(phi), cos(phi));
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
    // set state dimension
  n_x_ = 5;
  // set augmented dimension
  n_aug_ = n_x_ + 2;
  
  lambda_ = 3 - n_aug_;
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  // create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  //add measurement noise covariance matrix
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
       0                    , std_laspy_*std_laspy_;
  
  // add measurement noise covariance matrix
  R_radar_ = MatrixXd(3,3);
  R_radar_ <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
  is_initialized_ = false;
  //std::cout << "UKF instance initialized" << std::endl;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(false == is_initialized_)
  {
    
    float px;
    float py;
    float v;
    if(MeasurementPackage::RADAR ==  meas_package.sensor_type_)
    {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      px = rho*cos(phi);
      py = rho*sin(phi);
      v = rho_dot;
      
      P_ = MatrixXd::Identity(n_x_, n_x_);
      
      cout << "Initial RADAR measurement!" << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
      v = 0;
      
      P_ = MatrixXd::Identity(n_x_, n_x_);
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
      
      cout << "Initial LiDAR measurement!" << endl;
    }
    
    if(fabs(px) < 0.00005){
        px = 0.00005;
        cout << "intial x too small" << endl;
    }
    if(fabs(py) < 0.00005){
        py = 0.00005;
        cout << "intial y too small" << endl;
    }
    
    x_ << px, py, v, 0, 0;
    
    // set weights
    weights_.fill(0.5 / (n_aug_ + lambda_));
    weights_(0) = lambda_  / (lambda_ + n_aug_);
    
    is_initialized_ = true;
  }
  else
  {
    if(use_radar_ && MeasurementPackage::RADAR ==  meas_package.sensor_type_ )
    {
      R_ = R_radar_;
      this->Prediction((meas_package.timestamp_ - previous_timestamp_) / 1000000.0);
      this->UpdateRadar(meas_package);
    }
    else if(use_laser_ && MeasurementPackage::LASER ==  meas_package.sensor_type_ )
    {
      R_ = R_lidar_;
      this->Prediction((meas_package.timestamp_ - previous_timestamp_) / 1000000.0);
      this->UpdateLidar(meas_package);
    }
  }
  previous_timestamp_ = meas_package.timestamp_;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;
  
  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;
  
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  // print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  
  // predict sigma points
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    const double p_x = Xsig_aug(0,i);
    const double p_y = Xsig_aug(1,i);
    const double v = Xsig_aug(2,i);
    const double yaw = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  // print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
    
  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  
  // print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P_ << std::endl;
  //std::cout << "UKF::Prediction called" << std::endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  // LiDAR doesn't need linearization, update using linear Kalman equations
  // set measurement dimension, our LiDAR measures X and Y
  int n_z = 2;
  
  VectorXd z = VectorXd(n_z); 
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  
  MatrixXd H_ = MatrixXd(n_z, n_x_);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;
  
  //predicted measurement
  VectorXd z_pred = H_ * x_;
    
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // residual
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
  
  //std::cout << "UKF::UpdateLidar called" << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
    // set measurement dimension, radar can measure r, phi, and r_dot
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  /**
   * Student part begin
   */

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    const double p_x = Xsig_pred_(0,i);
    const double p_y = Xsig_pred_(1,i);
    const double v  = Xsig_pred_(2,i);
    const double yaw = Xsig_pred_(3,i);
    const double rho = sqrt((p_x*p_x) + (p_y*p_y));

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / max(0.000001,rho);  // r_dot
  }
  
  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred= Zsig* weights_;
  
//   mean predicted measurement
//   z_pred.fill(0.0);
//   for (int i=0; i < 2*n_aug_+1; ++i) {
//     z_pred = z_pred + weights_(i) * Zsig.col(i);
//   }
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    NormalizeAngle(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  S = S + R_;
  
    // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
    // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  //std::cout << "UKF::UpdateRadar called" << std::endl;

}
