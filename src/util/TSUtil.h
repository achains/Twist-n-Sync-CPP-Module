//
// Created by achains on 23.07.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_TSUTIL_H
#define TWIST_N_SYNC_CPP_MODULE_TSUTIL_H

#include "Eigen/Core"

#include <vector>

namespace TSUtil {
    // Implementation of numpy.arange function
    Eigen::VectorXd arangeEigen(double const & start, double const & finish, double const & step);

    // Implementation of numpy.diff function
    Eigen::VectorXd adjDiffEigen(Eigen::VectorXd const & data);

    // Implementation of scipy.signal.correlate function
    // arg1 : Eigen::VectorXd
    // arg2 : Eigen::VectorXd
    // return : Eigen::VectorXd
    Eigen::VectorXd eigenCrossCor(Eigen::VectorXd & data_1, Eigen::VectorXd & data_2);

    // Implementation of scipy.signal.correlate function
    // arg1 : std::vector<double>
    // arg2 : std::vector<double>
    // return : std::vector<double>
    std::vector<double> eigenCrossCor(std::vector<double> & data_1, std::vector<double> & data_2);

    Eigen::MatrixX3d vectorToEigMatrixX3d(std::vector<std::vector<double>> & data);

    Eigen::VectorXd vectorToEigVectorXd(std::vector<double> & data);

    Eigen::VectorXd getNormOfRows(Eigen::MatrixX3d const & data);

    struct CorrData{
        Eigen::VectorXd cross_cor;
        Eigen::Index initial_index;
    };
}


#endif //TWIST_N_SYNC_CPP_MODULE_TSUTIL_H
