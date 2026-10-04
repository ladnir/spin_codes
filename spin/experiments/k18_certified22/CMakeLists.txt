cmake_minimum_required(VERSION 3.20)
project(spin_k18_certified22 LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
get_filename_component(REPOSITORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(FROZEN "${REPOSITORY}/spin/experiments/k16_codesign_100us")
set(RETAINED "${REPOSITORY}/research/workstreams/k16_design/implementation")
add_executable(spin_k18_certified22 bench.cpp K18Setup.cpp K18Scalar.cpp K18Fast.cpp
  ${FROZEN}/ScalarAdapters.cpp
  ${FROZEN}/kernel/K16NativeSetup.cpp
  ${FROZEN}/outer_variants/OuterVariants.cpp
  ${RETAINED}/RsSetup.cpp ${RETAINED}/RsFast.cpp)
target_include_directories(spin_k18_certified22 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${FROZEN}" "${RETAINED}")
target_compile_features(spin_k18_certified22 PRIVATE cxx_std_20)
target_link_libraries(spin_k18_certified22 PRIVATE spin::spin)
target_compile_options(spin_k18_certified22 PRIVATE -mavx2 -mavx512f -mavx512vl
  -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=znver4 -fno-strict-aliasing -fno-lto)
enable_testing()
foreach(seed IN ITEMS 1 17 113)
  foreach(size IN ITEMS 2048 6144 65536 262144)
    add_test(NAME exact22_${size}_${seed} COMMAND spin_k18_certified22 ${size} ${seed} 0)
    set_tests_properties(exact22_${size}_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
  endforeach()
endforeach()
