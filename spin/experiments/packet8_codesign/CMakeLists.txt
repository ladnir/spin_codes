cmake_minimum_required(VERSION 3.20)
project(spin_packet8_codesign LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
get_filename_component(REPOSITORY "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(RETAINED "${REPOSITORY}/research/workstreams/k16_design/implementation")
add_executable(spin_packet8 bench.cpp Setup.cpp Scalar.cpp Fast.cpp OuterShared.cpp OuterHalf.cpp
  ${RETAINED}/RsSetup.cpp
  ../k16_codesign_100us/ScalarAdapters.cpp
  ../k16_codesign_100us/kernel/K16NativeSetup.cpp)
target_compile_features(spin_packet8 PRIVATE cxx_std_20)
target_link_libraries(spin_packet8 PRIVATE spin::spin)
target_include_directories(spin_packet8 PRIVATE "${RETAINED}")
if(MSVC)
  target_compile_options(spin_packet8 PRIVATE /arch:AVX512 /O2)
  target_compile_definitions(spin_packet8 PRIVATE SPIN_PACKET_ENABLE_VBMI=1)
else()
  target_compile_options(spin_packet8 PRIVATE -mavx2 -mavx512f -mavx512vl
    -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mprfchw -mtune=znver4
    -fno-strict-aliasing -fno-lto)
endif()
enable_testing()
foreach(size IN ITEMS 2048 6144 65536)
  foreach(seed IN ITEMS 1 17)
    add_test(NAME packet8_${size}_${seed} COMMAND spin_packet8 ${size} ${seed} 0)
    set_tests_properties(packet8_${size}_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
  endforeach()
endforeach()
