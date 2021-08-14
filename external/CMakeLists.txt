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

FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    set(eigen_EIGEN_BUILD_DOC OFF)
    set(eigen_BUILD_TESTS OFF)
    include_directories(twist_n_sync PUBLIC ${eigen_SOURCE_DIR})
endif()

FetchContent_GetProperties(spline)
if(NOT spline_POPULATED)
    FetchContent_Populate(spline)
    include_directories(twist_n_sync PUBLIC ${spline_SOURCE_DIR}/src)
endif()
