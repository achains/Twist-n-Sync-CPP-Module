include(FetchContent)

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
)

FetchContent_Declare(
        spline
        GIT_REPOSITORY https://github.com/ttk592/spline.git
        GIT_TAG origin/master
)

FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    set(eigen_EIGEN_BUILD_DOC OFF)
    set(eigen_BUILD_TESTS OFF)
    target_include_directories(twist_n_sync PUBLIC "${CMAKE_BINARY_DIR}/_deps/eigen-src")
endif()

FetchContent_GetProperties(spline)
if(NOT spline_POPULATED)
    FetchContent_Populate(spline)
    target_include_directories(twist_n_sync PUBLIC "${CMAKE_BINARY_DIR}/_deps/spline-src/src")
endif()
