#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//CONSTANT TURN RATE AND VELOCITY MODEL USED

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    is_initialized_ = false;
    previous_timestamp_ = 0;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    //Augmented mean vector
    x_aug = VectorXd(7);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Augmented Covariance Matrix
    P_aug = MatrixXd(7, 7);

    //Generated Sigma Points
    Xsig_aug_ = MatrixXd(7, 15); //MatrixXd(n_aug, 2*n_aug+1)

    //Predicted Sigma Points
    Xsig_pred_ = MatrixXd(5, 15); //MatrixXd(n_x, 2*n_aug+1)

/*********************
 * CHAGNGE THESESSESE*/
    // Process noise standard deviation longitudinal acceleration in m/s^2  CHANGE THESES
    std_a_ = 2.5;
/*********************
 * CHAGNGE THESESSESE*/
    // Process noise standard deviation yaw acceleration in rad/s^2   CHANGE THESES
    std_yawdd_ = 0.7;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15; //variance == std_laspx_^2


    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;//variance == std_laspy_^2

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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /***********************************************************
                      Initializtion Structure
    ************************************************************/
    if (!is_initialized_) {
        //Initialize x_, P_, previous_time, or any other variables
        x_.fill(0);
        n_x_ = 5;
        n_aug_ = 7;
        lambda_ = 3 - n_aug_;
        weights_ = VectorXd(2 * n_aug_ + 1);
        Xsig_aug_.fill(0);


        if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
             * Initialize State
             */
            //p_x position
            x_(0) = meas_package.raw_measurements_(0);
            //p_y position
            x_(1) = meas_package.raw_measurements_(1);
            //v
            x_(2) = 0;// actual value not known
            //psi
            x_(3) = 0;//actual value not known
            //psi_dot
            x_(4) = 0;//actual value not known
        } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            //initialize stuff here
            float rho = meas_package.raw_measurements_(0);
            float yaw = meas_package.raw_measurements_(1);
            float rho_dot = meas_package.raw_measurements_(2);
            x_(0) = rho * cos(yaw); //convert to px
            x_(1) = rho * sin(yaw); //convert to py
            x_(2) = rho_dot; //v approximation
            x_(3) = yaw;//yaw
            x_(4) = 0;//yaw dot
        }

        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }
    /***********************************************************
                      Control Structure
    ************************************************************/
    //Change in time
    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;

    Prediction(dt);

    //Update
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
    }
    //Reassigning Time
    previous_timestamp_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */
    /**************************************************************************
     *                  Step One Create Augmented Sigma Points
     *                  Lesson 7 Section 18
     * ************************************************************************/
    //set first five rows of augmented x equal to x
    cout << "Prediciton" << endl;
    x_aug.head(5) << x_;//make equal to x_

    x_aug(5) = 0.0;
    x_aug(6) = 0.0;


    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) << P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;


    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();


    //create augmented sigma points, number of sigma points= 2*n_aug+1 = 15

    //first column
    Xsig_aug_.col(0) = x_aug;

    for (int i = 0; i < n_aug_; i++) {
        //columns 2-8 or the first n_aug columns
        Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        //colums 9-15 or the second n_aug columns
        Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    /**************************************************************************
     *                  Predict Sigma Points
     *                  Lesson 7 Section 21
     * ************************************************************************/

    for (int i = 0; i < (2 * n_aug_ + 1); i++) {

        // extract values for better readability
        double p_x = Xsig_aug_(0, i);
        double p_y = Xsig_aug_(1, i);
        double v = Xsig_aug_(2, i);
        double yaw = Xsig_aug_(3, i);
        double yawd = Xsig_aug_(4, i);
        double nu_a = Xsig_aug_(5, i);
        double nu_yawdd = Xsig_aug_(6, i);

        //Passing sigma points through state prediction equations
        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        //predicted position
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        //predicted velocity and yaw/yaw rate
        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        //add noise
        //position noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);

        //velocity noise
        v_p = v_p + nu_a * delta_t;
        //yaw noise
        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        //Predicted sigma points, only five points because we do not return noise components
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;

    }



    /**************************************************************************
     *                  Predict Mean and Variance
     *                  Lesson 7 Section 24
     * ************************************************************************/

    // set weights
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int k = 1; k < (2 * n_aug_ + 1); k++) {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights_(k) = weight;
    }

    //predicted state mean
    x_.fill(0.0);
    for (int l = 0; l < (2 * n_aug_ + 1); l++) {
        x_ = x_ + weights_(l) * Xsig_pred_.col(l);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int m = 0; m < (2 * n_aug_ + 1); m++) {
        //state difference
        VectorXd x_diff = Xsig_pred_.col(m) - x_;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        P_ = P_ + weights_(m) * x_diff * x_diff.transpose();
    }

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
    cout << "Lidar Update" << endl;

    //    z(0) = x_(0);//px
    //    z(1) = x_(1);//py

    VectorXd z = VectorXd(2);
    z.fill(0);

    z<<meas_package.raw_measurements_;

    MatrixXd H_ = MatrixXd(2, 5);
    H_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0;
    MatrixXd R_ = MatrixXd(2, 2);
    R_ << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;

    VectorXd y = z - H_ * x_; //y = z-z_pred (z_pred = H_ * x_)
    MatrixXd Ht = H_.transpose();// 2x5 -> 5x2
    MatrixXd P_Ht = P_ * Ht;// 5x5 * 5x2 -> 5x2
    MatrixXd S = H_ * P_Ht+ R_; //2x5 * 5x2 + 2x2 =2x2
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;

    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    //new state
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;

    cout << "x" << endl;
    cout << x_ << endl;
    cout << endl;
    cout << "P" << endl;
    cout << P_ << endl;
    cout << endl;
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
    /**
     * Lesson 7 Section 27
     */
    cout << "Radar Update" << endl;

    MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);//matrix size 3x15
    //transform sigma points into measurement space
    for (int i = 0; i < (2 * n_aug_ + 1); i++) {//2n_aug+1 sigma points
        //extract values for better readability
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;//v in the x direction
        double v2 = sin(yaw) * v;//v in the y direction

        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                        //r
        Zsig(1, i) = atan2(p_y, p_x);                                 //phi
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(3);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) { //15 iterations 2*7+1
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S and cross correlation Matrix
    MatrixXd S = MatrixXd(3, 3);
    S.fill(0.0);
    MatrixXd Tc_ = MatrixXd(5, 3);
    Tc_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {//15 iterations
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }



    //add measurement noise covariance matrix
    MatrixXd R_ = MatrixXd(3, 3);
    R_.fill(0.0);
    R_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    S = S + R_;


    //Kalman gain K
    MatrixXd K_ = Tc_ * S.inverse();

    //Measurement Values
    VectorXd z_meas = VectorXd(3);
    z_meas = meas_package.raw_measurements_;

    //Residual
    VectorXd z_diff = z_meas - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    //update state mean and covariance matrix
    x_ = x_ + (K_ * z_diff);
    P_ = P_ - K_ * S * K_.transpose();

    std::cout << "Updated state x: " << std::endl << x_ << std::endl;
    std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}
