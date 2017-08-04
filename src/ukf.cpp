#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

//#define testingSig
//#define testingSigPred
//#define testingMeanPred
//#define testingRadarPred
//#define testingUKFUpdate

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
  for(int i=1; i<n_sig_; i++)
  {
    weights_(i) = 0.5/(lambda_+n_aug_);
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
  
#ifdef testingSig
  // Test Part
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;
  
  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //set example state
  VectorXd x = VectorXd(n_x);
  x <<   5.7441,
  1.3800,
  2.2049,
  0.5015,
  0.3528;
  
  //create example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
  -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
  0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
  -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
  -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
  
  // Change class variables to test variables
  x_ = x;
  P_ = P;
  std_a_ = std_a;
  std_yawdd_ = std_yawdd;
  lambda_ = lambda;
#endif
  
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

#ifdef testingSig
  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
#endif
  
  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd sigmaPoints, double delta_t)
{
#ifdef testingSigPred
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug <<
  5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
  1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
  2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
  0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
  0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
  0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
  0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

  delta_t = 0.1; //time diff in sec
  n_x_ = n_x;
  n_aug_ = n_aug;
  sigmaPoints = Xsig_aug;
#endif
  
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

#ifdef testingSigPred
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred = Xsig_pred_;
  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
#endif

  return;
}

void UKF::PredictMeanAndCovariance()
{
#ifdef testingMeanPred
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  n_x_ = n_x;
  n_aug_ = n_aug;
  lambda_ = lambda;
  Xsig_pred_ = Xsig_pred;
  
#endif
  
  x_.fill(0.0);
  for(int i=0; i<n_sig_; ++i)
  {
    x_ += Xsig_pred_.col(i)*weights_(i);
  }
  
  P_.fill(0.0);
  for(int i=0; i<n_sig_; ++i)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI)
    {
      x_diff(3) -= 2.*M_PI;
    }
    while (x_diff(3)<-M_PI)
    {
      x_diff(3) += 2.*M_PI;
    }
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  
#ifdef testingMeanPred
  //create vector for predicted state
  VectorXd x = VectorXd(n_x);
  
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);
  
  x = x_;
  P = P_;
  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;
#endif
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
  z_pred.fill(0.0);
  
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
  S.fill(0.0);
  
  for (int i=0; i<n_sig_; ++i)
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI)
    {
      z_diff(1)-=2.*M_PI;
    }
    while (z_diff(1)< -M_PI)
    {
      z_diff(1)+=2.*M_PI;
    }
    
    S = S + weights_(i) * z_diff * z_diff.transpose();
    //S += weights_(i)*(Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S += R;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
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
  
#ifdef testingRadarPred
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  //int n_z = 3;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //radar measurement noise standard deviation radius in m
  double std_radr = 0.3;
  
  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;
  
  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;
  
  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  n_x_ = n_x;
  n_aug_ = n_aug;
  lambda_ = lambda;
  std_radr_ = std_radr;
  std_radphi_ = std_radphi;
  std_radrd_ = std_radrd;
  Xsig_pred_ = Xsig_pred;
//  weights_ = weights;
  
#endif
  
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
  z_pred.fill(0.0);
  
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
  S.fill(0.0);
  
  for (int i=0; i<n_sig_; ++i)
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI)
    {
     z_diff(1)-=2.*M_PI;
    }
    while (z_diff(1)< -M_PI)
    {
     z_diff(1)+=2.*M_PI;
    }
    
    S = S + weights_(i) * z_diff * z_diff.transpose();
    
    //S += weights_(i)*(Zsig.col(i) - z_pred)*(Zsig.col(i) - z_pred).transpose();
  }
  S += R;
  
#ifdef testingRadarPred
  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;
#endif
  
#ifdef testingUKFUpdate
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z = 3;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
  5.93637,
  1.49035,
  2.20528,
  0.536853,
  0.353577;
  
  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
  -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
  -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;
  
  //create example matrix with sigma points in measurement space
  Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
  6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
  0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
  2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;
  
  //create example vector for mean predicted measurement
  z_pred = VectorXd(n_z);
  z_pred <<
  6.12155,
  0.245993,
  2.10313;
  
  //create example matrix for predicted measurement covariance
  S = MatrixXd(n_z,n_z);
  S <<
  0.0946171, -0.000139448,   0.00407016,
  -0.000139448,  0.000617548, -0.000770652,
  0.00407016, -0.000770652,    0.0180917;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
  5.9214,   //rho in m
  0.2187,   //phi in rad
  2.0062;   //rho_dot in m/s
  
  n_x_ = n_x;
  n_aug_ = n_aug;
  lambda_ = lambda;
  Xsig_pred_ = Xsig_pred;
  x_ = x;
  P_ = P;
  meas_package.raw_measurements_ = z;
  
#endif
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for(int i=0; i<n_sig_; ++i)
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI)
    {
      z_diff(1) -= 2.*M_PI;
    }
    while (z_diff(1)<-M_PI)
    {
      z_diff(1)+=2.*M_PI;
    }
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI)
    {
      x_diff(3)-=2.*M_PI;
    }
    while (x_diff(3)<-M_PI)
    {
      x_diff(3)+=2.*M_PI;
    }
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    //Tc += weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
  P_ = P_-K*S*K.transpose();
  
#ifdef testingUKFUpdate
  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
#endif
  
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
