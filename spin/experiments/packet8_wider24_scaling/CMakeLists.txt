cmake_minimum_required(VERSION 3.20)
project(spin_packet8_wider24_scaling LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
set(base ../packet8_wider24)
add_executable(spin_packet8_wider24_scaling bench.cpp StoreRoute.cpp StoreOuter.cpp
  ../packet8_wider24_opt/RandomizerVariants.cpp
  ${base}/Setup.cpp ${base}/Scalar.cpp ${base}/Fast.cpp ${base}/Outer.cpp)
target_compile_features(spin_packet8_wider24_scaling PRIVATE cxx_std_20)
target_link_libraries(spin_packet8_wider24_scaling PRIVATE spin::spin)
if(MSVC)
  target_compile_options(spin_packet8_wider24_scaling PRIVATE /arch:AVX512 /O2)
  target_compile_definitions(spin_packet8_wider24_scaling PRIVATE SPIN_PACKET_ENABLE_VBMI=1)
else()
  target_compile_options(spin_packet8_wider24_scaling PRIVATE -mavx2 -mavx512f
    -mavx512vl -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mprfchw
    -mtune=znver4 -fno-strict-aliasing -fno-lto)
endif()
enable_testing()
foreach(size IN ITEMS 2048 6144 65536 262144 1048576)
  foreach(seed IN ITEMS 1 17)
    add_test(NAME packet8_scaling_${size}_${seed}
      COMMAND spin_packet8_wider24_scaling ${size} ${seed} 0 99 normal 0)
    set_tests_properties(packet8_scaling_${size}_${seed}
      PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
  endforeach()
endforeach()
