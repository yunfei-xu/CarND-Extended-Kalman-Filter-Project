#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /**
   * TODO:
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  // Process
  ekf_.F_ = Eigen::MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
    0, 1, 0, 1,
    0, 0, 1, 0,
    0, 0, 0, 1;

  ekf_.P_ = Eigen::MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;

  noise_ax_ = 9;
  noise_ay_ = 9;

  // Measurement
  H_laser_ = Eigen::MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;

  R_laser_ = Eigen::MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
    0, 0.0225;

//  Hj_ = Eigen::MatrixXd(3, 4);

  R_radar_ = Eigen::MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF(){ }

void FusionEKF::processMeasurement(const MeasurementPackage& measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
     * TODO:
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * Remember: you'll need to convert radar from polar to cartesian coordinates.
     */
    // first measurement
    std::cout << "EKF: " << std::endl;
    ekf_.x_ = Eigen::VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordinates and initialize state.
       */
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho * std::cos(phi), rho * std::sin(phi), 0, 0;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       * Initialize state.
       */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    std::cout << "EKF initialized!" << std::endl;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   * TODO:
   * Update the state transition matrix F according to the new elapsed time.
   *  - Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2)       = dt;
  ekf_.F_(1, 3)       = dt;

  // set the process covariance matrix Q
  float dt2 = dt * dt;
  float dt3 = dt * dt2;
  float dt4 = dt * dt3;
  ekf_.Q_ = Eigen::MatrixXd(4, 4);
  ekf_.Q_ << dt4 / 4 * noise_ax_, 0, dt3 / 2 * noise_ax_, 0,
    0, dt4 / 4 * noise_ay_, 0, dt3 / 2 * noise_ay_,
    dt3 / 2 * noise_ax_, 0, dt2 * noise_ax_, 0,
    0, dt3 / 2 * noise_ay_, 0, dt2 * noise_ay_;

  ekf_.predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   * TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools_.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.updateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
} // FusionEKF::processMeasurement
