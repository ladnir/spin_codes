cmake_minimum_required(VERSION 3.20)
project(spin_packet_k16 LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
set(SPIN_TUNE znver4 CACHE STRING "Compiler scheduling target")
add_executable(spin_packet_k16 bench.cpp CachedRoute.cpp CachedOuter.cpp)
target_link_libraries(spin_packet_k16 PRIVATE spin::spin)
if(MSVC)
  target_compile_options(spin_packet_k16 PRIVATE /arch:AVX512)
  target_compile_definitions(spin_packet_k16 PRIVATE SPIN_PACKET_ENABLE_VBMI=1)
else()
  target_compile_options(spin_packet_k16 PRIVATE -mavx2 -mavx512f -mavx512vl
    -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=${SPIN_TUNE}
    -fno-strict-aliasing -fno-lto)
endif()
