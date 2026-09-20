cmake_minimum_required(VERSION 3.20)
project(SPINCodes VERSION 0.1.0 LANGUAGES CXX)

option(SPIN_BUILD_TESTS "Build standalone API tests" ON)
if(SPIN_BUILD_TESTS)
    enable_testing()
endif()
add_subdirectory(spin)
