enable_testing()

# ###############################################
# Fetching googletest
include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG origin/main
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)
# ###############################################

set(TEST_SOURCES final_result_tests.cpp final_result_tests.h main.cpp tsutil_tests.cpp tsutil_tests.h)

# Copying data from inputData to binary dir
file(GLOB SENSORS_DATA ${CMAKE_SOURCE_DIR}/tests/inputData/sensors)
file(COPY ${SENSORS_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(run_tests ${TEST_SOURCES})

# Including fetched content
target_include_directories(run_tests PUBLIC "${CMAKE_BINARY_DIR}/_deps/googletest-src/googletest/include")

target_link_libraries(run_tests PUBLIC gtest_main)
target_link_libraries(run_tests PUBLIC twist_n_sync)

include(GoogleTest)
gtest_discover_tests(run_tests)