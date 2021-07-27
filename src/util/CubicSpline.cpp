//
// Created by achains on 26.07.2021.
//

#include "CubicSpline.h"
#include "TSUtil.h"

#include <iostream>

#include <memory>

CubicSpline::CubicSpline(Eigen::VectorXd const & X, Eigen::VectorXd const & Y):
    spline_(std::vector<double>(X.data(), X.data() + X.size()),
            std::vector<double>(Y.data(), Y.data() + Y.size())){}

CubicSpline::CubicSpline(std::vector<double> const & X, std::vector<double> const & Y):
    spline_(X, Y, tk::spline::cspline_hermite) {}

void CubicSpline::calculateDerivative() {
    derivative_.resize(spline_.get_x().size());
    size_t i = 0;
    for (auto& elem: spline_.get_x()){
        derivative_[i++] = spline_.deriv(1, elem);
    }

    derivative_.front()= derivative_.back() = 0.0;
}

Eigen::Matrix4Xd CubicSpline::getCoefficients() {
    if (derivative_.empty()){
        calculateDerivative();
    }
    std::vector<double> Y = spline_.get_y();
    Eigen::VectorXd A = TSUtil::vectorToEigVectorXd(Y);
    Eigen::VectorXd B = TSUtil::vectorToEigVectorXd(derivative_);

    assert(A.size() == B.size());

    double dx = spline_.get_x()[1] - spline_.get_x()[0];

    Eigen::VectorXd C(A.size() - 1);
    Eigen::VectorXd D(A.size() - 1);

    for (Eigen::Index i = 0; i < A.size() - 1; ++i){
        C[i] = -(2.0 * B[i] + B[i + 1]) / dx + 3.0 * (A[i + 1] - A[i]) / (dx * dx);
        D[i] = -2.0 * C[i] / (3. * dx) + (B[i + 1] - B[i]) / (3 * dx * dx);
    }

    Eigen::Matrix4Xd coefficients(4, A.size() - 1);

    // There are N - 1 segments, where N is number of knots.
    // C, D already have A.size() - 1, where A.size() equals N
    // Getting rid of the last redundant element in A and B.
    coefficients << A.transpose().block(0, 0, A.size() - 1, 1),
                    B.transpose().block(0, 0, A.size() - 1, 1),
                    C.transpose(), D.transpose();

    return coefficients;
}