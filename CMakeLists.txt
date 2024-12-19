cmake_minimum_required(VERSION 3.28)
project(heart_beat_generation)

set(CMAKE_CXX_STANDARD 17)

find_package(Python3 COMPONENTS Development REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(heart_beat_generation src/mcmc_wrapper.c src/likelihood.c src/likelihood.h src/util.c src/util.h)
target_include_directories(heart_beat_generation PRIVATE ${Python3_INCLUDE_DIRS})
target_link_libraries(heart_beat_generation PUBLIC OpenMP::OpenMP_CXX)
