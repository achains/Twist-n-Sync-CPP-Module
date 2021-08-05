set(SOURCES TimeSync.cpp TimeSync.h)

add_library(twist_n_sync SHARED ${SOURCES})

target_include_directories(twist_n_sync PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
