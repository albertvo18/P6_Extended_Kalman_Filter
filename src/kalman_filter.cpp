#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Constructor
KalmanFilter::KalmanFilter() {}

// Destructor
KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  // Predict
  // Predicted (apriori) state estimate
  // Predicted (apriori) estimate covariance
  // x_  : object state    
  // P_  : object covariance matrix    
  // F_  : state transition matrix    


  // H_  : measurement matrix    
  // R_  : measurement covariance matrix    
  // Q_  : process covariance matrix    
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

// Kalman Filter predict function
void KalmanFilter::Predict() {
  x_ = F_ * x_ ;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // calculate error and update 
  VectorXd y = z - H_ * x_;
  Update_with_Measurement(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);


  

  double rho = sqrt(px*px + py*py);
  // check if rho is 0 to avoid dividing by 0
  // if rho is 0, then return without updating
  if ( rho == 0)
    return;

  double theta = atan2(py, px);
  double rho_dot = (px*vx + py*vy) / rho;

  // initialize h  
  VectorXd h = VectorXd(3);
  h << rho, theta, rho_dot;

  // Calculate y (Measurement) 
  VectorXd y = z - h;
  while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    } else {
      y(1) += M_PI;
    }
  }
  Update_with_Measurement(y);
}

void KalmanFilter::Update_with_Measurement(const VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
