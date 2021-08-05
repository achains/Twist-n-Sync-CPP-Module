cmake_minimum_required(VERSION 3.14)
project(Twist_n_Sync_CPP_Module)

# Set this flag ON if you want to build tests
option(BUILD_TESTS "Build an executable with tests" OFF)

# Further key-section adds build with google sanitizers
option(USE_ASAN "Build with ASan sanitizer" OFF)

# Compiler prerequisites
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# By default select Debug build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

# If Build type is Release, then turn on optimizations
if (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

if (USE_ASAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
endif()

# Prevent In-Source Builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
message(FATAL_ERROR "In-source build detected!")
endif()

add_subdirectory(src)
add_subdirectory(src/util)
add_subdirectory(external)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()