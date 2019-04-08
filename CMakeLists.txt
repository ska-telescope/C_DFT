cmake_minimum_required(VERSION 3.10.2)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message(">>> Building project in DEBUG mode")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

IF(CMAKE_BUILD_TYPE MATCHES Release)
    message(">>> Building project in RELEASE mode")
ENDIF(CMAKE_BUILD_TYPE MATCHES Release)

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Base direct fourier transform project
project(dft)
add_executable(dft main.cpp direct_fourier_transform.c)
target_link_libraries(dft m)

# Unit testing for dft
project(tests)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(tests unit_testing.cpp direct_fourier_transform.c)
target_link_libraries(tests ${GTEST_LIBRARIES} pthread)
