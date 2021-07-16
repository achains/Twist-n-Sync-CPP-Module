//
// Created by achains on 18.07.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H
#define TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H

#include "Eigen/Core"
#include <vector>

class TimeSync {
 public:
    explicit TimeSync(Eigen::MatrixX3d  gyro_first, Eigen::MatrixX3d  gyro_second,
                      Eigen::VectorXd  ts_first, Eigen::VectorXd  ts_second, bool do_resample);
 private:
    // 3D angular velocities from devices' gyros
    Eigen::MatrixX3d gyro_first_;
    Eigen::MatrixX3d gyro_second_;

    // Gyros' timestamps
    Eigen::VectorXd ts_first_;
    Eigen::VectorXd ts_second_;

    // Flag to do resampling of angular velocities
    bool do_resample_;

    // Note: np.diff consumes the first element, but the std::adj_diff -- not!
    static Eigen::ArrayXd adjDiffEigen(Eigen::ArrayXd const & A);

    static Eigen::VectorXd arangeEigen(double const & start, double const & finish, double const & step);

    // Took this function from https://forum.kde.org/viewtopic.php?f=74&t=118619
    // TODO: Overview this function and fix warnings
    static std::vector<double> eigenCrossCor(std::vector<double>& xCorrInputVecFirst,
                                             std::vector<double>& xCorrInputVecSecond);

    void obtainDelay();

    void resample(double const & accuracy);
};


#endif //TWIST_N_SYNC_CPP_MODULE_TIMESYNC_H
