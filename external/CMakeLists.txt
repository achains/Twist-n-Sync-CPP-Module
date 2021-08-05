include(FetchContent)

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

FetchContent_MakeAvailable(eigen)
FetchContent_MakeAvailable(spline)

# Including fetched content
target_include_directories(twist_n_sync PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
target_include_directories(twist_n_sync PUBLIC "${CMAKE_BINARY_DIR}/_deps/spline-src/src")
