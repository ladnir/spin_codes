cmake_minimum_required(VERSION 3.20)
project(spin_locality LANGUAGES CXX)
set(SPIN_TUNE znver4 CACHE STRING "Scheduling target")
get_filename_component(SPIN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../spin" ABSOLUTE)
add_subdirectory("${SPIN_ROOT}" spin)
add_executable(locality experiment.cpp)
target_link_libraries(locality PRIVATE spin::spin)
target_include_directories(locality PRIVATE "${SPIN_ROOT}/src/kernels")
target_compile_definitions(locality PRIVATE SPIN_BCH_AVX512=1 SPIN_WORKSPACE_ROUTING_OPT=1)
# Mechanical load-layout specialization of the unchanged four-row circuit.
# Keep all arithmetic, sharing, and fixed-width SIMD operations intact.
file(READ "${SPIN_ROOT}/src/kernels/generated/BchAvx512.cpp" MAPPED_BCH)
string(REGEX REPLACE "_mm512_loadu_si512\\(a\\+4\\*([0-9]+)\\)" "_mm512_loadu_si512(a+permutation[\\1])" MAPPED_BCH "${MAPPED_BCH}")
string(REPLACE "part0(const block* a,block* x)" "part0(const block* a,block* x,const unsigned* permutation)" MAPPED_BCH "${MAPPED_BCH}")
string(REPLACE "bchTranspose4(const block* a,block* x)" "bchTranspose4Mapped(const block* a,block* x,const unsigned* permutation)" MAPPED_BCH "${MAPPED_BCH}")
string(REPLACE "part0(a,x);" "part0(a,x,permutation);" MAPPED_BCH "${MAPPED_BCH}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/MappedBch.cpp" "${MAPPED_BCH}")
add_library(locality_bch OBJECT "${CMAKE_CURRENT_BINARY_DIR}/MappedBch.cpp")
foreach(STRIDE IN ITEMS 16 20)
  file(READ "${SPIN_ROOT}/src/kernels/generated/BchAvx512.cpp" STRIDED_BCH)
  string(REPLACE "_mm512_loadu_si512(a+4*" "_mm512_loadu_si512(a+${STRIDE}*" STRIDED_BCH "${STRIDED_BCH}")
  string(REPLACE "bchTranspose4(" "bchTranspose4Stride${STRIDE}(" STRIDED_BCH "${STRIDED_BCH}")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/BchStride${STRIDE}.cpp" "${STRIDED_BCH}")
  target_sources(locality_bch PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/BchStride${STRIDE}.cpp")
endforeach()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/schedule_bch.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/resynth_bch.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/gfni_bch.py"
  "${SPIN_ROOT}/src/kernels/generated/BchCircuit.h"
  "${SPIN_ROOT}/src/kernels/generated/BchAvx512.cpp")
foreach(SCHEDULE IN ITEMS Greedy Release Reverse Restrict)
  set(SCHEDULE_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/Bch${SCHEDULE}.cpp")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/schedule_bch.py"
    "${SPIN_ROOT}/src/kernels/generated/BchAvx512.cpp"
    "${SPIN_ROOT}/src/kernels/generated/BchCircuit.h" "${SCHEDULE}"
    OUTPUT_FILE "${SCHEDULE_SOURCE}" COMMAND_ERROR_IS_FATAL ANY)
  target_sources(locality_bch PRIVATE "${SCHEDULE_SOURCE}")
endforeach()
foreach(SHARING IN ITEMS 6 8)
  set(SHARING_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/BchShare${SHARING}.cpp")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/resynth_bch.py"
    "${SPIN_ROOT}/src/kernels/generated/BchCircuit.h"
    "${SPIN_ROOT}/../research/scripts/probe_bch_forward_xor_circuit.py" "${SHARING}"
    OUTPUT_FILE "${SHARING_SOURCE}" COMMAND_ERROR_IS_FATAL ANY)
  target_sources(locality_bch PRIVATE "${SHARING_SOURCE}")
endforeach()
target_compile_features(locality_bch PRIVATE cxx_std_20)
set(GFNI_SOURCE "${CMAKE_CURRENT_BINARY_DIR}/BchGfni.cpp")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/gfni_bch.py"
  "${SPIN_ROOT}/src/kernels/generated/BchCircuit.h"
  OUTPUT_FILE "${GFNI_SOURCE}" COMMAND_ERROR_IS_FATAL ANY)
target_sources(locality_bch PRIVATE "${GFNI_SOURCE}")
if(NOT MSVC)
  set_source_files_properties("${GFNI_SOURCE}" PROPERTIES COMPILE_OPTIONS "-mgfni;-mavx512bw")
endif()
target_include_directories(locality_bch PRIVATE "${SPIN_ROOT}/src/kernels" "${SPIN_ROOT}/include")
target_sources(locality PRIVATE $<TARGET_OBJECTS:locality_bch>)
if(MSVC)
  target_compile_options(locality PRIVATE /arch:AVX512)
  target_compile_options(locality_bch PRIVATE /arch:AVX512)
else()
  # Match Fast.cpp, including its larger register file for the 19-word state.
  # This isolated executable already requires AVX-512; public dispatch is unchanged.
  target_compile_options(locality PRIVATE -mavx2 -mavx512f -mavx512vl -mtune=${SPIN_TUNE} -fno-strict-aliasing)
  target_compile_options(locality_bch PRIVATE -mavx2 -mavx512f -mavx512vl -mtune=${SPIN_TUNE} -fno-strict-aliasing)
endif()
