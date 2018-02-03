#include "kalman_filter.h"
#include <math.h>
#include <iostream>

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter(){
  I_ = Eigen::MatrixXd::Identity(4, 4);
}

KalmanFilter::~KalmanFilter(){ }

void KalmanFilter::init(Eigen::VectorXd& x, Eigen::MatrixXd& P,
                        Eigen::MatrixXd& F, Eigen::MatrixXd& Q,
                        Eigen::MatrixXd& H, Eigen::MatrixXd& R)
{
  x_ = x;
  P_ = P;
  F_ = F;
  Q_ = Q;
  H_ = H;
  R_ = R;
}

void KalmanFilter::predict()
{
  /**
   * TODO:
   * predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::update(const Eigen::VectorXd& z)
{
  /**
   * TODO:
   * update the state by using Kalman Filter equations
   */
  Eigen::VectorXd y = z - H_ * x_;
  Eigen::MatrixXd S = H_ * P_ * H_.transpose() + R_;
  Eigen::MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ += K * y;
  P_  = (I_ - K * H_) * P_;
}

void KalmanFilter::updateEKF(const Eigen::VectorXd& z)
{
  /**
   * TODO:
   * update the state by using Extended Kalman Filter equations
   */

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  Eigen::VectorXd hx(3);
  hx(0) = std::sqrt(px * px + py * py);
  hx(1) = std::atan2(py, px);
  hx(2) = (px * vx + py * vy) / hx(0);

  Eigen::VectorXd y = z - hx;
  if (y(1) > M_PI) {
    y(1) -= 2 * M_PI;
  } else if (y(1) < -M_PI) {
    y(1) += 2 * M_PI;
  }
  Eigen::MatrixXd S = H_ * P_ * H_.transpose() + R_;
  Eigen::MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ += K * y;
  P_  = (I_ - K * H_) * P_;
}
