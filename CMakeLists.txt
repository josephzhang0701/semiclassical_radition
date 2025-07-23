cmake_minimum_required(VERSION 3.11)
project(semiclassical_radiation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)

# 包含 GSL 头文件，链接 GSL 库目录
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

# 收集当前目录下所有 .cpp 源文件
file(GLOB COMMON_SRCS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

# 排除两个 main 源文件，留作各自的入口
list(REMOVE_ITEM COMMON_SRCS
        "${CMAKE_CURRENT_SOURCE_DIR}/main_r.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/multi_particle_tracer.cpp"
)

# -------------------------------------------------------------------
# 旧可执行文件 irun：使用 main_r.cpp
add_executable(irun
        main_r.cpp
        ${COMMON_SRCS}
)
target_link_libraries(irun
        PRIVATE ${GSL_LIBRARIES} m
)

# -------------------------------------------------------------------
# 新可执行文件 mpt：使用 multi_particle_tracer.cpp
add_executable(mpt
        multi_particle_tracer.cpp
        ${COMMON_SRCS}
)
target_link_libraries(mpt
        PRIVATE ${GSL_LIBRARIES} m
)