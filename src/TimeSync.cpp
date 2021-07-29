//
// Created by achains on 18.07.2021.
//

#include "TimeSync.h"
#include "util/TSUtil.h"
#include "util/CubicSpline.h"

#include "Eigen/Dense"

#include <numeric>

TimeSync::TimeSync(std::vector<std::vector<double>> & gyro_first,
                   std::vector<std::vector<double>> & gyro_second,
                   std::vector<double> & ts_first,
                   std::vector<double> & ts_second,
                   bool const & do_resample):
                        gyro_first_(std::move(TSUtil::vectorToEigMatrixX3d(gyro_first))),
                        gyro_second_(std::move(TSUtil::vectorToEigMatrixX3d(gyro_second))),
                        ts_first_(std::move(TSUtil::vectorToEigVectorXd(ts_first))),
                        ts_second_(std::move(TSUtil::vectorToEigVectorXd(ts_second))),
                        do_resample_(do_resample) {}


TSUtil::CorrData TimeSync::getInitialIndex() const {

    Eigen::VectorXd norm_first = TSUtil::getNormOfRows(gyro_first_);
    Eigen::VectorXd norm_second = TSUtil::getNormOfRows(gyro_second_);

    Eigen::VectorXd cross_cor = TSUtil::eigenCrossCor(norm_first, norm_second);

    return {cross_cor, std::distance(cross_cor.begin(), std::max(cross_cor.begin(), cross_cor.end()))};
}

Eigen::MatrixX3d TimeSync::interpolateGyro(Eigen::VectorXd const & ts_old, Eigen::MatrixX3d const & gyro_old,
                                           Eigen::VectorXd const & ts_new) {
    Eigen::MatrixX3d gyro_new(gyro_old.rows(), 3);
    gyro_new << TSUtil::interpolate(ts_old, gyro_old(Eigen::all, 0), ts_new),
                TSUtil::interpolate(ts_old, gyro_old(Eigen::all, 1), ts_new),
                TSUtil::interpolate(ts_old, gyro_old(Eigen::all, 2), ts_new);
    return gyro_new;
}

void TimeSync::resample(double const & accuracy){
    double time_first_mean = TSUtil::adjDiffEigen(ts_first_).mean();
    double time_second_mean = TSUtil::adjDiffEigen(ts_second_).mean();

    double dt = std::min({accuracy, time_first_mean, time_second_mean});

    if (do_resample_){
        Eigen::VectorXd ts_first_new = TSUtil::arangeEigen(ts_first_[0], *ts_first_.end() + dt, dt);
        Eigen::VectorXd ts_second_new = TSUtil::arangeEigen(ts_second_[0], *ts_second_.end() + dt, dt);

        gyro_first_ = TimeSync::interpolateGyro(ts_first_, gyro_first_, ts_first_new);
        gyro_second_ = TimeSync::interpolateGyro(ts_second_, gyro_second_, ts_second_new);
    }
}

void TimeSync::obtainDelay(){
    // Correction of index numbering
    Eigen::Index shift = -gyro_first_.rows() + 1;

    // Cross-cor estimation
    TSUtil::CorrData corr_data = TimeSync::getInitialIndex();
    corr_data.initial_index += shift;

    Eigen::MatrixX3d tmp_xx1 = gyro_first_;
    Eigen::MatrixX3d tmp_xx2 = gyro_second_;

    if (corr_data.initial_index > 0){
        tmp_xx1 = gyro_first_(Eigen::seq(0, Eigen::last - corr_data.initial_index), Eigen::all);
        tmp_xx2 = gyro_second_(Eigen::seq(corr_data.initial_index, Eigen::last), Eigen::all);
    }
    else if (corr_data.initial_index < 0){
        tmp_xx1 = gyro_first_(Eigen::seq(-corr_data.initial_index, Eigen::last), Eigen::all);
        tmp_xx2 = gyro_second_(Eigen::seq(0, corr_data.initial_index), Eigen::all);
    }

    size_t size = std::min(gyro_first_.rows(), gyro_second_.rows());
    tmp_xx1 = tmp_xx1(Eigen::seq(0, size - 1), Eigen::all);
    tmp_xx2 = tmp_xx2(Eigen::seq(0, size - 1), Eigen::all);

    // Calibration
    Eigen::MatrixX3d M = (tmp_xx2.transpose() * tmp_xx1) * (tmp_xx1.transpose() * tmp_xx1).inverse();
    gyro_first_ = (M * gyro_first_.transpose()).transpose();

    // Cross-correlation re-estimation
    corr_data = TimeSync::getInitialIndex();

    // Cross-cor, based cubic spline coefficients
    CubicSpline cubic_spline(TSUtil::arangeEigen(0., static_cast<double>(corr_data.cross_cor.size())),
                                corr_data.cross_cor);

    Eigen::Matrix4Xd spline_coefficients = cubic_spline.getCoefficients();

    Eigen::VectorXd coeffs = spline_coefficients(Eigen::all, corr_data.initial_index);

    // Check cubic spline derivative sign and redefine initial_index if needed
    if (coeffs(Eigen::last - 1) < 0) {
        corr_data.initial_index -= 1;
        coeffs = spline_coefficients(Eigen::all, corr_data.initial_index);
    }

    // Solve quadratic equation to obtain roots
    Eigen::Index order = coeffs.size() - 1;
    Eigen::VectorXd equation(order);
    for (auto i = 0; i < order; ++i){
        equation[i] = static_cast<double>(order - i) * coeffs[i];
    }
    Eigen::VectorXd roots = TSUtil::quadraticRoots(equation);


    auto result = *std::max(roots.begin(), roots.end());
    std::vector<double> check_solution(order);
    for (int i = 0; i < order; ++i)
        check_solution[i] = static_cast<double>(order - i) * coeffs[i] *
                std::pow((roots[0] + roots[1]) / 2, (order - i - 1));
    if (std::accumulate(check_solution.begin(), check_solution.end(), 0.0) < 0.0)
        result = *std::min(roots.begin(), roots.end());

    time_delay_ = result;
}

double TimeSync::getTimeDelay() const {
    return time_delay_;
}
