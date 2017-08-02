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
  // For intiliasing the filter
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  // Matrix sizes
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_+1;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  lambda_ = 3 - n_aug_;
  
  Xsig_pred_ = MatrixXd(n_x_,n_sig_);
  
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i=1; i<n_sig_; ++i)
  {
    weights_(i) = 1/(2*(lambda_+n_aug_));
  }
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  // INITILISATION
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Conver Rader from polar to cartesian coordinates and initilise state
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float px = rho*cos(theta);
      float py = rho*sin(theta);
      x_ << px,py,0,0,0;
    }
    else
    {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;
    }
    
    // Initialise P
    P_ <<  1,0,0,0,0,
           0,1,0,0,0,
           0,0,1,0,0,
           0,0,0,1,0,
           0,0,0,0,1;
    // TIME??????
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  else if (use_laser_)
  {
    UpdateLidar(meas_package);
  }
  
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd sigmaPoints = GenerateSigmaPoints();
  PredictSigmaPoints(sigmaPoints, delta_t);
  PredictMeanAndCovariance();
}

MatrixXd UKF::GenerateSigmaPoints()
{
  // Create augmented mean
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  // Create augmented covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  // Create square root
  MatrixXd sqrtP = P_aug.llt().matrixL();
  
  // Create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_; ++i)
  {
    Xsig_aug.col(i+1)         = x_aug + sqrt(lambda_+n_aug_) * sqrtP.col(i);
    Xsig_aug.col(i+1+n_aug_)  = x_aug - sqrt(lambda_+n_aug_) * sqrtP.col(i);
  }
  
  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd sigmaPoints, double delta_t)
{
  for(int i=0; i<n_sig_; ++i)
  {
    double px = sigmaPoints(0,i);
    double py = sigmaPoints(1,i);
    double v = sigmaPoints(2,i);
    double yaw = sigmaPoints(3,i);
    double yawd = sigmaPoints(4,i);
    double nu_a = sigmaPoints(n_x_,i);
    double nu_yawdd = sigmaPoints(n_x_+1,i);
    
    VectorXd x = VectorXd(n_x_);
    x << px,py,v,yaw,yawd;
    
    VectorXd change = VectorXd(n_x_);
    if (fabs(yawd) > 0.001)
    {
      change << v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw)),
                v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw)),
                0, yawd*delta_t,0;
    }
    else
    {
      change << v*cos(yaw)*delta_t,
                v*sin(yaw)*delta_t,
                0,0,0;
    }
    
    VectorXd noise = VectorXd(n_x_);
    noise << 0.5*delta_t*delta_t*cos(yaw)*nu_a,
            0.5*delta_t*delta_t*sin(yaw)*nu_a,
            delta_t*nu_a,
            0.5*delta_t*delta_t*nu_yawdd,
            delta_t*nu_yawdd;
    
    Xsig_pred_.col(i) = x+change+noise;
  }
  return;
}

void UKF::PredictMeanAndCovariance()
{
  for(int i=0; i<n_sig_; ++i)
  {
    x_ += Xsig_pred_.col(i)*weights_(i);
    
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while(x_diff(3)>M_PI)
    {
      x_diff(3) -= 2.*M_PI;
    }
    while(x_diff(3) <-M_PI)
    {
      x_diff(3) += 2.*M_PI;
    }
    
    P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
  }
  return;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Measurement dimension
  int n_z = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z,n_sig_);
  
  //transform sigma points into measurement space
  for (int i=0; i<n_sig_; ++i)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    
    Zsig.col(i) << px, py;
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //calculate mean predicted measurement
  for (int i=0; i<n_sig_; ++i)
  {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  
  R <<  std_laspx_*std_laspx_,0,
        0,std_laspy_*std_laspy_;
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for (int i=0; i<n_sig_; ++i)
  {
    S += weights_(i)*(Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S += R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  for(int i=0; i<n_sig_; ++i)
  {
    Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
  P_ = P_-K*S*K.transpose();
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  // Measurement dimension
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z,n_sig_);
  
  //transform sigma points into measurement space
  for (int i=0; i<n_sig_; ++i)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);
    
    double rho = sqrt(px*px + py*py);
    double phi = atan2(py,px);
    double rhod = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;
    
    Zsig.col(i) << rho, phi, rhod;
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //calculate mean predicted measurement
  for (int i=0; i<n_sig_; ++i)
  {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  
  R <<  std_radr_*std_radr_,0,0,
        0,std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for (int i=0; i<n_sig_; ++i)
  {
    S += weights_(i)*(Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S += R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  for(int i=0; i<n_sig_; ++i)
  {
    Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
  P_ = P_-K*S*K.transpose();
  
  return;
}

//void UKF::UKFUpdate(MatrixXd S, MatrixXd Zsig, MatrixXd z_pred, MeasurementPackage meas_package)
//{
//  //create matrix for cross correlation Tc
//  MatrixXd Tc = MatrixXd(n_x_, n_z);
//  
//  //calculate cross correlation matrix
//  for(int i=0; i<n_sig_; ++i)
//  {
//    Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
//  }
//  
//  //calculate Kalman gain K;
//  MatrixXd K = Tc * S.inverse();
//  
//  //update state mean and covariance matrix
//  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
//  P_ = P_-K*S*K.transpose();
//  
//  return;
//}
