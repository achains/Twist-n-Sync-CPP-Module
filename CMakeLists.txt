cmake_minimum_required(VERSION 3.19)
project(Twist_n_Sync_CPP_Module)

# By default select Debug build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# compiler settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(lib)

set(SOURCE_FILES main.cpp src/TimeSync.cpp)
set(EXECUTABLE_NAME app_debug)

add_executable(${EXECUTABLE_NAME} ${SOURCE_FILES})

target_link_libraries(${EXECUTABLE_NAME} PRIVATE alglib)
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")