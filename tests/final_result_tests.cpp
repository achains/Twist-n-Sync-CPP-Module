//
// Created by achains on 04.08.2021.
//

#include "final_result_tests.h"
#include "gtest/gtest.h"
#include "tsutil_tests.h"
#include "TimeSync.h"

#include <fstream>


namespace final_result_tests{
    void csvParser(std::vector<std::vector<double>> & gyro, std::vector<double> & time, char const * fileName){
        std::ifstream openFile(fileName);
        EXPECT_TRUE(openFile);

        std::string line{};
        while (std::getline(openFile, line)){
            std::vector<double> gyro_dimension(3);
            for (int i = 0; i < 3; ++i){
                gyro_dimension[i] = std::stod(line.substr(0, line.find(',')));
                line = line.substr(line.find(',') + 1);
            }
            gyro.push_back(gyro_dimension);
            time.push_back(std::stod(line));
        }
        openFile.close();
    }

    TEST(ParserTest, Random5Values){
        std::vector<std::vector<double>> gyro;
        std::vector<double> time;

        csvParser(gyro, time, "sensors/gyro_0-session_0.csv");
        // Lines 40-44
        std::vector<std::vector<double>> expected_gyro =
                {{0.9608788,-2.5956225,-0.033207856}, {0.94255286,-2.6231115,0.007720115},
                 {0.8936837,-2.6921394,0.08041308,6}, {0.86558384,-2.7226827,0.10973461},
                 {0.81427115,-2.762389,0.15677124}};
        std::vector<double> expected_time = {
                690644511255257, 690644513254990, 690644515254723,690644517254523, 690644519254323
        };
        std::vector<std::vector<double>> sliced_gyro(gyro.begin() + 39, gyro.begin() + 44);
        std::vector<double> sliced_time(time.begin() + 39, time.begin() + 44);

        EXPECT_TRUE(tsutil_tests::compareVectors(sliced_time, expected_time, 1e-8));
        for (int i = 39; i < 44; ++i)
            EXPECT_TRUE(tsutil_tests::compareVectors(sliced_gyro[i], expected_gyro[i], 1e-8));
    }

    TEST(TimeSyncTest, TimeSyncResult){
        std::vector<std::vector<double>> gyro_first, gyro_second;
        std::vector<double> time_first, time_second;
        csvParser(gyro_first, time_first, "sensors/gyro_0-session_2.csv");
        csvParser(gyro_second, time_second, "sensors/gyro_1-session_2.csv");

        for (double& elem : time_first) elem /= 1e9;
        for (double& elem: time_second) elem /= 1e9;

        TimeSync timesync(gyro_first, gyro_second, time_first, time_second, false);
        timesync.resample(1.0);
        timesync.obtainDelay();

//        std::clog.precision(17);
//        std::clog << "Delay: " << timesync.getTimeDelay() << std::endl;
    }

}