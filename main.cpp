#include "TimeSync.h"
#include "util/CubicSpline.h"

#include <iostream>
#include <fstream>

#include <vector>
#include <string>

using std::vector;
using std::string;

void collect_data(vector<vector<double>> & gyro, vector<double> & time, char const * fileName){
    std::ifstream file(fileName);

    string line;
    while (std::getline(file, line)){
        vector<double> gyro_dimension(3);
        for (int i = 0; i < 3; ++i){
            gyro_dimension[i] = std::stod(line.substr(0, line.find(',')));
            line = line.substr(line.find(',') + 1);
        }
        gyro.push_back(gyro_dimension);
        time.push_back(std::stod(line));
    }
    std::clog << ">>> " << fileName << " successfully parsed" << std::endl;
    file.close();
}

int main(){

    std::cout.precision(17);

    vector<vector<double>> gyro_1;
    vector<vector<double>> gyro_2;

    vector<double> time_1;
    vector<double> time_2;

    char const * FIRST_GYRO_DATA = "sensors/gyro_0-session_1.csv";
    char const * SECOND_GYRO_DATA = "sensors/gyro_1-session_1.csv";

    collect_data(gyro_1, time_1, FIRST_GYRO_DATA);
    collect_data(gyro_2, time_2, SECOND_GYRO_DATA);

    std::clog << ">>> gyro_1.size() = " << gyro_1.size() << " gyro_2.size() = " << gyro_2.size() << std::endl;

    for (auto& elem: time_1) elem /= 1e9;
    for (auto& elem: time_2) elem /= 1e9;

    TimeSync time_sync(gyro_1, gyro_2, time_1, time_2, false);

    time_sync.resample(1.0);
    time_sync.obtainDelay();

    std::clog << ">>> time_delay = " << time_sync.getTimeDelay() << std::endl;

    return 0;
}