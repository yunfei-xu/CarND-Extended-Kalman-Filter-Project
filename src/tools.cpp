#include <iostream>
#include "tools.h"

Tools::Tools(){ }

Tools::~Tools(){ }

Eigen::VectorXd Tools::CalculateRMSE(const std::vector<Eigen::VectorXd>& estimations,
                                     const std::vector<Eigen::VectorXd>& ground_truth)
{
  Eigen::VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  size_t n = estimations.size();

  // Check the validity of the following inputs:
  // * the estimation vector size should not be zero
  // * the estimation vector size should equal ground truth vector size
  assert(n != 0);
  assert(n == ground_truth.size());

  // Accumulate squared residuals
  for (size_t i = 0; i < n; i++) {
    Eigen::VectorXd delta = estimations[i] - ground_truth[i];
    Eigen::VectorXd delta2 = delta.array() * delta.array();
    rmse = rmse + delta2;
  }

  // Calculate the mean
  rmse /= n;

  // Calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

Eigen::MatrixXd Tools::CalculateJacobian(const Eigen::VectorXd& x_state)
{
  Eigen::MatrixXd Hj(3, 4);

  // Recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float dist2 = px * px + py * py;
  float dist  = std::sqrt(dist2);
  float dist3 = dist2 * dist;

  // Check division by zero
  assert(dist2 >= 1e-6);

  // Compute the Jacobian matrix
  Hj << px / dist, py / dist, 0, 0,
    -py / dist2, px / dist2, 0, 0,
    py * (vx * py - vy * px) / dist3, px * (vy * px - vx * py) / dist3, px / dist, py / dist;

  return Hj;
}
