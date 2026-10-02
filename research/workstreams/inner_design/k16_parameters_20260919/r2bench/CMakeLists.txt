cmake_minimum_required(VERSION 3.20)
project(SpinK16TwoRoundExperiment LANGUAGES CXX)
# This overlay is deliberately isolated from the reusable encoder defaults.
set(SPIN_BIDIRECTIONAL_SOURCE "" CACHE PATH "" FORCE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../../spin_optimized" base)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/overlay.py"
  "${CMAKE_CURRENT_BINARY_DIR}/base/generated" COMMAND_ERROR_IS_FATAL ANY)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS overlay.py)
set(CMAKE_CXX_STANDARD 20)
add_executable(r2_recurrence recurrence.cpp)
target_link_libraries(r2_recurrence PRIVATE spin_half_transpose)
enable_testing()
add_test(NAME two_round_adjoint COMMAND r2_recurrence)
