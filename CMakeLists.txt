cmake_minimum_required(VERSION 3.17)
project(Twist_n_Sync_CPP_Module)

# By default select Debug build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# compiler settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add_compile_options("$<$<CONFIG:Debug>:-O2;-DELPP_THREAD_SAFE;-g>")

set(SOURCE_FILES main.cpp src/TimeSync.cpp src/util/TSUtil.cpp src/util/CubicSpline.cpp)
set(EXECUTABLE_NAME app_debug)

# Adding .csv gyro data files
file(GLOB SENSORS_DATA ${CMAKE_SOURCE_DIR}/tests/inputData/sensors)
file(COPY ${SENSORS_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})

target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib/spline/src")
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/")