cmake_minimum_required(VERSION 3.20)
project(SpinOptimized LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC OR NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
  message(FATAL_ERROR "This integration requires x86-64 GCC or Clang")
endif()
option(SPIN_BCH_AVX512 "Enable runtime-dispatched four-row BCH" OFF)
option(SPIN_TEST_NO_AVX512 "Test-only: suppress AVX-512 detection" OFF)
option(SPIN_BUILD_BENCHMARK "Build the serial, locked experiment harness" OFF)
option(SPIN_BUILD_WIDE "Import Hypercat's 256/512-bit forward kernels" OFF)
option(SPIN_FORWARD_FOUR "Experimental, unvalidated four-row forward circuit" OFF)
option(SPIN_GENERAL_LENGTHS "Naturally aligned lengths for the transposed library" ON)
set(SPIN_DIRECT_MAX_K 458752 CACHE STRING "Upper K for runtime-length AVX-512 direct routing; 0 disables")
option(SPIN_FORWARD_PROFILE "Diagnostic forward-stage clocks" OFF)
option(SPIN_STAGE_PROFILE "Diagnostic tiled-kernel clocks and Linux perf counters" OFF)
set(SPIN_BCH_SCHEDULE global CACHE STRING "Four-row BCH schedule: global, local8, local16, local32, single, reload")
if(NOT SPIN_BCH_SCHEDULE MATCHES "^(global|local8|local16|local32|single|reload)$")
  message(FATAL_ERROR "Unknown SPIN_BCH_SCHEDULE")
endif()
set(SPIN_BIDIRECTIONAL_SOURCE "" CACHE PATH "Optional upstream native/spin source directory (copied, never edited)")
set(SPIN_TUNE generic CACHE STRING "GCC/Clang tuning target; changes scheduling, not required instructions")
option(SPIN_K16_DIRECT "Use direct routing for K=2^16 with the AVX-512 backend" ON)
option(SPIN_K18_DIRECT "Use direct routing for K=2^18 with the AVX-512 backend" ON)
set(SPIN_K20_EXPERIMENT 0 CACHE STRING "K20 research override: 0=tiled default, 1=direct32, 2=direct24, 3/4=direct prefetch16/64, 5=stream stores, 6/7=batched stream/cached stores")
if(NOT SPIN_K20_EXPERIMENT MATCHES "^[0-7]$")
  message(FATAL_ERROR "SPIN_K20_EXPERIMENT must be an integer from 0 through 7")
endif()
set(SPIN_K18_EXPERIMENT auto CACHE STRING "K18 routing: auto=normal, 0=tiled, 1=direct32, 2=direct24, 3/4=prefetch16/64")
set(K18_MODE "${SPIN_K18_EXPERIMENT}")
if(K18_MODE STREQUAL "auto")
  set(K18_MODE 0)
  if(SPIN_K18_DIRECT)
    set(K18_MODE 1)
  endif()
endif()
set(SPIN_SMALL_EXPERIMENT auto CACHE STRING "Research override: auto=normal, 0=control, 1=32-bit, 2=packed24, 3/4=prefetch16/64")
set(SMALL_MODE "${SPIN_SMALL_EXPERIMENT}")
if(SMALL_MODE STREQUAL "auto")
  set(SMALL_MODE 0)
  if(SPIN_K16_DIRECT)
    set(SMALL_MODE 1)
  endif()
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
mark_as_advanced(SPIN_TEST_NO_AVX512 SPIN_GENERAL_LENGTHS SPIN_DIRECT_MAX_K
  SPIN_K16_DIRECT SPIN_K18_DIRECT SPIN_SMALL_EXPERIMENT SPIN_K18_EXPERIMENT
  SPIN_K20_EXPERIMENT SPIN_BCH_SCHEDULE SPIN_STAGE_PROFILE
  SPIN_BUILD_WIDE SPIN_BIDIRECTIONAL_SOURCE SPIN_FORWARD_FOUR SPIN_FORWARD_PROFILE)
get_filename_component(ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(GEN "${CMAKE_CURRENT_BINARY_DIR}/generated")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/generate.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/bch_local.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/k16.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/k16_r2.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/k18.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/k20_tiled.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/stage_profile.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/StageProfile.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/ImtRounds.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/K16R2Inner.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/k16_map.py"
  "${CMAKE_CURRENT_SOURCE_DIR}/K16Inner.h"
  "${ROOT}/workstreams/inner_design/bch_20260919/generate.py"
  "${ROOT}/scripts/probe_bch_forward_xor_circuit.py"
  "${ROOT}/workstreams/bare_bch_rm2sub/Spin.h"
  "${ROOT}/workstreams/bare_bch_rm2sub/Inner.h"
  "${ROOT}/workstreams/bare_bch_rm2sub/WorkspaceRouting.h"
  "${ROOT}/workstreams/bare_bch_rm2sub/generated/SelectedMaps.h"
  "${ROOT}/workstreams/bare_bch_rm2sub/generated/BchCircuit.h"
  "${ROOT}/workstreams/bare_bch_rm2sub/generated/BchCircuit.cpp"
  "${ROOT}/workstreams/inner_design/asymmetric/bch256/weight5/implementation/Weight5Spin.cpp"
  "${ROOT}/workstreams/inner_design/asymmetric/bch256/weight5/implementation/Weight5Inner.h"
  "${ROOT}/workstreams/inner_design/asymmetric/bch256/weight5/implementation/AsymmetricMap.h"
  "${ROOT}/workstreams/rate_quarter_bch/implementation/generated/QuarterCircuit.h")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generate.py" "${GEN}"
  COMMAND_ERROR_IS_FATAL ANY)
if(NOT SPIN_BCH_SCHEDULE STREQUAL "global")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/bch_local.py" "${GEN}" "${SPIN_BCH_SCHEDULE}" COMMAND_ERROR_IS_FATAL ANY)
endif()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/small.py")
if(NOT SMALL_MODE STREQUAL "0")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/small.py" "${GEN}" "${SMALL_MODE}" COMMAND_ERROR_IS_FATAL ANY)
endif()
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/k16.py" "${GEN}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/k16_r2.py" "${GEN}" COMMAND_ERROR_IS_FATAL ANY)
if(NOT K18_MODE STREQUAL "0")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/k18.py" "${GEN}" "${K18_MODE}" COMMAND_ERROR_IS_FATAL ANY)
endif()
set(VENDOR "${ROOT}/constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation/vendor")
if(SPIN_K20_EXPERIMENT MATCHES "^[567]$")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/k20_tiled.py" "${GEN}" "${SPIN_K20_EXPERIMENT}" COMMAND_ERROR_IS_FATAL ANY)
elseif(NOT SPIN_K20_EXPERIMENT STREQUAL "0")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/k18.py" "${GEN}" "${SPIN_K20_EXPERIMENT}" 20 COMMAND_ERROR_IS_FATAL ANY)
endif()
# Keep all ordinary code at the old AVX2 baseline, even on an AVX-512 build host.
if(SPIN_STAGE_PROFILE)
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/stage_profile.py" "${GEN}" COMMAND_ERROR_IS_FATAL ANY)
endif()
# Do not propagate AVX-512 options to callers or allow LTO to cross the ISA boundary.
add_library(spin_baseline INTERFACE)
target_compile_options(spin_baseline INTERFACE -march=x86-64-v3 "-mtune=${SPIN_TUNE}" -mpclmul -mvpclmulqdq -fno-lto)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lengths.py")
if(SPIN_GENERAL_LENGTHS)
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/lengths.py" "${GEN}" "${SPIN_DIRECT_MAX_K}" COMMAND_ERROR_IS_FATAL ANY)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/generic.py" "${CMAKE_CURRENT_SOURCE_DIR}/GenericSpin.h")
  execute_process(COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/generic.py" "${GEN}" COMMAND_ERROR_IS_FATAL ANY)
endif()
add_library(spin_half_transpose STATIC "${GEN}/Spin.cpp" "${GEN}/generated/BchCircuit.cpp")
if(SPIN_GENERAL_LENGTHS)
  target_sources(spin_half_transpose PRIVATE "${GEN}/GenericSpin.cpp")
endif()
target_include_directories(spin_half_transpose PUBLIC "${GEN}" "${VENDOR}")
target_link_libraries(spin_half_transpose PUBLIC spin_baseline)
target_compile_definitions(spin_half_transpose PUBLIC SPIN_GENERAL_LENGTHS=$<BOOL:${SPIN_GENERAL_LENGTHS}>)
target_compile_definitions(spin_half_transpose PRIVATE
  SPIN_BCH_AVX512=$<BOOL:${SPIN_BCH_AVX512}> SPIN_TEST_NO_AVX512=$<BOOL:${SPIN_TEST_NO_AVX512}>
  SPIN_GROUPED_A=1 SPIN_WORKSPACE_ROUTING_OPT=1 W5_MASKED=0 W5_SHARED_EMISSION=0 W5_WRITE_PREFETCH=0 W5_PREFETCH=32)
if(SPIN_BCH_AVX512)
  add_library(spin_bch512 OBJECT "${GEN}/generated/BchAvx512.cpp" "${GEN}/Fast.cpp")
  target_include_directories(spin_bch512 PRIVATE "${GEN}" "${VENDOR}")
  target_compile_options(spin_bch512 PRIVATE -march=x86-64-v3 "-mtune=${SPIN_TUNE}" -mpclmul -mvpclmulqdq -mavx512f -mavx512vl -fno-lto)
  target_compile_definitions(spin_bch512 PRIVATE SPIN_GROUPED_A=1 W5_MASKED=0 W5_SHARED_EMISSION=0 W5_WRITE_PREFETCH=0 W5_PREFETCH=32)
  set_property(TARGET spin_bch512 PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
  target_sources(spin_half_transpose PRIVATE $<TARGET_OBJECTS:spin_bch512>)
endif()
set_property(TARGET spin_half_transpose PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
include(CTest)
if(BUILD_TESTING)
  if(SPIN_GENERAL_LENGTHS)
    add_executable(spin_generic_test generic_test.cpp)
    target_link_libraries(spin_generic_test PRIVATE spin_half_transpose)
    add_test(NAME generic_elements COMMAND spin_generic_test)
    set_tests_properties(generic_elements PROPERTIES RUN_SERIAL TRUE)
    add_executable(spin_lengths_test lengths_test.cpp)
    target_link_libraries(spin_lengths_test PRIVATE spin_half_transpose)
    add_test(NAME aligned_lengths COMMAND spin_lengths_test)
    set_tests_properties(aligned_lengths PROPERTIES RUN_SERIAL TRUE)
  endif()
  if(SPIN_BCH_AVX512)
    add_executable(spin_bch_basis bch_basis_test.cpp)
    target_link_libraries(spin_bch_basis PRIVATE spin_half_transpose)
    add_test(NAME bch_basis COMMAND spin_bch_basis)
    set_tests_properties(bch_basis PROPERTIES SKIP_RETURN_CODE 77 RUN_SERIAL TRUE)
  endif()
  add_test(NAME k16_map_exact COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_k16_map.py")
  add_executable(spin_correctness "${ROOT}/workstreams/inner_design/asymmetric/bch256/weight5/implementation/correctness.cpp")
  target_link_libraries(spin_correctness PRIVATE spin_half_transpose)
  add_test(NAME dense_reference COMMAND spin_correctness)
  add_executable(spin_backends backend_test.cpp)
  target_link_libraries(spin_backends PRIVATE spin_half_transpose)
  add_test(NAME backend_equivalence COMMAND spin_backends)
  add_executable(spin_k16_test k16_test.cpp)
  target_link_libraries(spin_k16_test PRIVATE spin_half_transpose)
  add_test(NAME k16_transpose COMMAND spin_k16_test)
  set_tests_properties(k16_transpose PROPERTIES RUN_SERIAL TRUE)
  set_tests_properties(dense_reference backend_equivalence PROPERTIES RUN_SERIAL TRUE)
endif()
if(SPIN_BUILD_BENCHMARK)
  add_executable(spin_benchmark benchmark.cpp)
  target_link_libraries(spin_benchmark PRIVATE spin_half_transpose)
endif()
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/OptionalForward.cmake")
