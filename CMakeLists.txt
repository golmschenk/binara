cmake_minimum_required(VERSION 3.26.5)
project(binara)

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 COMPONENTS Development REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(binara src/mcmc_wrapper.c src/likelihood.c src/likelihood.h src/util.c src/util.h)
target_include_directories(binara PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(binara PUBLIC OpenMP::OpenMP_CXX m)
