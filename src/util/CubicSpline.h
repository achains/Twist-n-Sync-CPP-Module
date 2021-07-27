//
// Created by achains on 26.07.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_CUBICSPLINE_H
#define TWIST_N_SYNC_CPP_MODULE_CUBICSPLINE_H

#include "spline.h"
#include "Eigen/Core"

class CubicSpline {
 public:
    explicit CubicSpline(Eigen::VectorXd const & X, Eigen::VectorXd const & Y);
    explicit CubicSpline(std::vector<double> const & X, std::vector<double> const & Y);

    Eigen::Matrix4Xd getCoefficients();

 private:
    void calculateDerivative();

    tk::spline spline_;
    std::vector<double> derivative_;
};


#endif //TWIST_N_SYNC_CPP_MODULE_CUBICSPLINE_H
