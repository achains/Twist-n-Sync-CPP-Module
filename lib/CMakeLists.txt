include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG origin/master
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG origin/master
)

FetchContent_Declare(
        spline
        GIT_REPOSITORY https://github.com/ttk592/spline.git
        GIT_TAG origin/master
)


FetchContent_MakeAvailable(googletest)
FetchContent_MakeAvailable(eigen)
FetchContent_MakeAvailable(spline)

# Including fetched content
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_BINARY_DIR}/_deps/googletest-src/googletest/include")
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
target_include_directories(${EXECUTABLE_NAME} PUBLIC "${CMAKE_BINARY_DIR}/_deps/spline-src/src")
