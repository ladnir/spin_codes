cmake_minimum_required(VERSION 3.20)
project(spin_k16_codesign LANGUAGES CXX)

# Isolated research build: consume an existing library snapshot rather than
# rebuilding a concurrently edited production tree.
find_package(spin CONFIG REQUIRED)
get_filename_component(REPOSITORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(RETAINED "${REPOSITORY}/research/workstreams/k16_design/implementation")
# Preserve the retained sources. Only adapt the retired public comparison enum
# and its label in a generated copy of the research driver.
file(READ "${RETAINED}/driver.cpp" DRIVER)
string(REPLACE "PacketT64S16" "PacketRsT64S20" DRIVER "${DRIVER}")
string(REPLACE "bch" "production" DRIVER "${DRIVER}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/retained_driver.cpp" "${DRIVER}")
add_executable(spin_k16_rs
  ${CMAKE_CURRENT_BINARY_DIR}/retained_driver.cpp ${RETAINED}/RsSetup.cpp ${RETAINED}/RsScalar.cpp
  ${RETAINED}/RsFast.cpp ${RETAINED}/RsLarge.cpp)
target_include_directories(spin_k16_rs PRIVATE "${RETAINED}")
target_compile_features(spin_k16_rs PRIVATE cxx_std_20)
target_link_libraries(spin_k16_rs PRIVATE spin::spin)
target_compile_options(spin_k16_rs PRIVATE -mavx2 -mavx512f -mavx512vl
  -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=znver4
  -fno-strict-aliasing -fno-lto)
enable_testing()
add_test(NAME retained_rs_correctness COMMAND spin_k16_rs check)
set_tests_properties(retained_rs_correctness PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_k16_codesign bench.cpp ScalarAdapters.cpp kernel/K16Setup.cpp
  kernel/K16Kernels.cpp kernel/K16T128.cpp kernel/K16NativeSetup.cpp
  kernel/K16NativeOuter.cpp kernel/K16Monomial.cpp kernel/K16Paired.cpp
  kernel/K16PairedOptimized.cpp kernel/K16PairedOptimizedSetup.cpp
  kernel/K16PairedWide.cpp
  kernel/K16PairedBasis.cpp kernel/K16PairedBasisSetup.cpp
  kernel/K16PairedBasisRaw.cpp
  kernel/K16Paired15.cpp kernel/K16Paired15Setup.cpp kernel/K16Paired15Scalar.cpp
  kernel/K16PairedShuffle.cpp kernel/K16PairedShuffleSetup.cpp
  kernel/K16Paired15Shuffle.cpp kernel/K16Paired15ShuffleSetup.cpp
  kernel/K16PairedNoAlias.cpp
  kernel/K16Paired15NoAlias.cpp
  kernel/K16PairedFold.cpp kernel/K16Paired15Fold.cpp
  outer_variants/OuterVariants.cpp
  outer_variants/OuterShuffleVariants.cpp
  outer_variants/OuterNonTemporal.cpp
  ${RETAINED}/RsSetup.cpp ${RETAINED}/RsFast.cpp)
target_compile_features(spin_k16_codesign PRIVATE cxx_std_20)
target_link_libraries(spin_k16_codesign PRIVATE spin::spin)
target_include_directories(spin_k16_codesign PRIVATE "${RETAINED}")
target_compile_options(spin_k16_codesign PRIVATE -mavx2 -mavx512f -mavx512vl
  -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=znver4
  -fno-strict-aliasing -fno-lto)
foreach(mode RANGE 0 52)
  foreach(size IN ITEMS 2048 6144 65536)
    foreach(seed IN ITEMS 1 17)
      add_test(NAME candidate_${mode}_${size}_${seed} COMMAND spin_k16_codesign ${mode} ${size} ${seed} 0)
      set_tests_properties(candidate_${mode}_${size}_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
    endforeach()
  endforeach()
endforeach()
