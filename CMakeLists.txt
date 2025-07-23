cmake_minimum_required(VERSION 3.11)
project(semiclassical_radiation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)

# 包含 GSL 头文件，链接 GSL 库目录
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

# -------------------------------------------------------------------
# 收集除了两个 main 以外的所有 .cpp
file(GLOB_RECURSE SRC_COMMON CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)
list(REMOVE_ITEM SRC_COMMON
        "${CMAKE_CURRENT_SOURCE_DIR}/main_r.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/multi_particle_tracer.cpp"
)

# -------------------------------------------------------------------
# 定义一个小宏，自动添加可执行 target + run target
macro(add_app NAME MAIN_SRC)
    # 可执行文件
    add_executable(${NAME}
            ${MAIN_SRC}
            ${SRC_COMMON}
    )
    target_link_libraries(${NAME} PRIVATE ${GSL_LIBRARIES} m)

    # 一键运行 target
    add_custom_target(run_${NAME}
            COMMAND ${NAME}
            DEPENDS ${NAME}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            COMMENT "Run ${NAME}"
    )
endmacro()

# -------------------------------------------------------------------
# 旧版入口：main_r.cpp → irun
add_app(irun main_r.cpp)

# 新版多粒子追踪：multi_particle_tracer.cpp → mpt
add_app(mpt multi_particle_tracer.cpp)