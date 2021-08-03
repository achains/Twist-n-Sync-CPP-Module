cmake_minimum_required(VERSION 3.17)
project(Twist_n_Sync_CPP_Module)

# By default select Debug build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# compiler settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_FILES src/TimeSync.cpp src/util/TSUtil.cpp src/util/CubicSpline.cpp tests/tests.cpp tests/tsutil_tests.cpp)
set(EXECUTABLE_NAME app_debug)

# Adding .csv gyro data files
file(GLOB SENSORS_DATA ${CMAKE_SOURCE_DIR}/tests/inputData/sensors)
file(COPY ${SENSORS_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})

add_subdirectory("src")
add_subdirectory("lib")

target_link_libraries(${EXECUTABLE_NAME} PUBLIC gtest_main)