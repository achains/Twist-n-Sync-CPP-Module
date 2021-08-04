//
// Created by achains on 04.08.2021.
//

#ifndef TWIST_N_SYNC_CPP_MODULE_FINAL_RESULT_TESTS_H
#define TWIST_N_SYNC_CPP_MODULE_FINAL_RESULT_TESTS_H

#include "Eigen/Core"

#include <vector>

namespace final_result_tests{
    // Stores parsed from fileName gyro data to arg1 and time data to arg2
    void csvParser(std::vector<std::vector<double>> & gyro, std::vector<double> & time, char const * fileName);
}

#endif //TWIST_N_SYNC_CPP_MODULE_FINAL_RESULT_TESTS_H