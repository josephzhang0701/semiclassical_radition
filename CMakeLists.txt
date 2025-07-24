cmake_minimum_required(VERSION 3.11)
project(semiclassical_radiation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# GSL via pkg-config
find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)

# OpenMP 支持（使用 Homebrew 安装的 clang/clang++）
find_package(OpenMP REQUIRED)

# 包含目录
include_directories(
        ${GSL_INCLUDE_DIRS}
)
link_directories(${GSL_LIBRARY_DIRS})

# 收集源文件
file(GLOB COMMON_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
list(REMOVE_ITEM COMMON_SRCS
        "${CMAKE_CURRENT_SOURCE_DIR}/main_r.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/multi_particle_tracer.cpp"
)

# irun
add_executable(irun
        main_r.cpp
        ${COMMON_SRCS}
)
target_link_libraries(irun
        PRIVATE
        ${GSL_LIBRARIES}
        m
        OpenMP::OpenMP_CXX
)

# mpt
add_executable(mpt
        multi_particle_tracer.cpp
        ${COMMON_SRCS}
)
target_link_libraries(mpt
        PRIVATE
        ${GSL_LIBRARIES}
        m
        OpenMP::OpenMP_CXX
)