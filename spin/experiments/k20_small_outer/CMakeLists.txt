cmake_minimum_required(VERSION 3.20)
project(spin_k20_small_outer LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
get_filename_component(REPOSITORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(EXACT "${REPOSITORY}/spin/experiments/k18_certified22")
set(FROZEN "${REPOSITORY}/spin/experiments/k16_codesign_100us")
set(RETAINED "${REPOSITORY}/research/workstreams/k16_design/implementation")
add_library(exact22_objects OBJECT
  ${EXACT}/K18Setup.cpp ${EXACT}/K18Scalar.cpp ${EXACT}/K18Fast.cpp
  ${FROZEN}/ScalarAdapters.cpp ${FROZEN}/kernel/K16NativeSetup.cpp
  ${FROZEN}/outer_variants/OuterVariants.cpp
  ${RETAINED}/RsSetup.cpp ${RETAINED}/RsFast.cpp
  K20Route.cpp K20Outer.cpp)
target_include_directories(exact22_objects PUBLIC "${EXACT}" "${FROZEN}" "${RETAINED}")
target_compile_features(exact22_objects PUBLIC cxx_std_20)
target_link_libraries(exact22_objects PUBLIC spin::spin)
target_compile_options(exact22_objects PUBLIC -mavx2 -mavx512f -mavx512vl
  -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=znver4 -fno-strict-aliasing -fno-lto)
add_executable(spin_k20_small_outer bench.cpp)
target_link_libraries(spin_k20_small_outer PRIVATE exact22_objects)
# Reuse the unchanged literal-map checker, including routed-state comparison.
add_executable(spin_k20_exact22_checks ${EXACT}/bench.cpp)
target_link_libraries(spin_k20_exact22_checks PRIVATE exact22_objects)
add_executable(spin_k20_store_checks check.cpp)
target_link_libraries(spin_k20_store_checks PRIVATE exact22_objects)
enable_testing()
foreach(seed IN ITEMS 1 41 113)
  add_test(NAME exact22_k20_${seed} COMMAND spin_k20_exact22_checks 1048576 ${seed} 0)
  set_tests_properties(exact22_k20_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
  add_test(NAME stores_k20_${seed} COMMAND spin_k20_store_checks 1048576 ${seed})
  set_tests_properties(stores_k20_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
endforeach()
foreach(k IN ITEMS 2048 65536 262144 264192 1046528 1050624)
  add_test(NAME stores_k${k} COMMAND spin_k20_store_checks ${k} 41)
  set_tests_properties(stores_k${k} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
endforeach()
