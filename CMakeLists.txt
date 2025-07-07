cmake_minimum_required(VERSION 3.11)
project(semiclassical_radition)

set(CMAKE_CXX_STANDARD 11)
find_package(PkgConfig REQUIRED)
pkg_search_module(GSL REQUIRED gsl)

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

file(GLOB SRC "*.cpp")
add_executable(irun.x ${SRC})
target_link_libraries(irun.x ${GSL_LIBRARIES} m)