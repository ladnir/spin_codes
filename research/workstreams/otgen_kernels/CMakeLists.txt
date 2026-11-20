cmake_minimum_required(VERSION 3.20)
project(otgen_kernels CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SPIN_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../spin" CACHE PATH "SPIN source tree")
set(OTGEN_ARCH "-march=native" CACHE STRING "ISA flags for the experiment kernels")
option(OTGEN_FUSED "Build the fused SPIN+AES experiment (needs an AVX-512 VAES host)" ON)
separate_arguments(OTGEN_ARCH_LIST UNIX_COMMAND "${OTGEN_ARCH}")

set(SPIN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(SPIN_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
add_subdirectory(${SPIN_SOURCE_DIR} spin EXCLUDE_FROM_ALL)

add_executable(otgen_bench bench.cpp)
target_link_libraries(otgen_bench PRIVATE spin::spin)
target_compile_options(otgen_bench PRIVATE ${OTGEN_ARCH_LIST} -O3)

if(OTGEN_FUSED)
  # Patched copies of SPIN's two AVX-512 transpose stages with an inlined hook.
  # The library sources are only read; the copies live in the build tree.
  include(${CMAKE_CURRENT_SOURCE_DIR}/patch_spin.cmake)
  otgen_patch(${SPIN_SOURCE_DIR}/src/certified22/RouteFast.cpp ${CMAKE_BINARY_DIR}/hooked/RouteFast.cpp
    "for(std::size_t epoch=first;--epoch;) {"
    "for(std::size_t epoch=first;--epoch;) { OTGEN_HOOK_ROUTE();")
  otgen_patch(${SPIN_SOURCE_DIR}/src/rs8/Fast.cpp ${CMAKE_BINARY_DIR}/hooked/Fast.cpp
    "group<Stream>(in+groupStride*g,out+128*g,plan.field.data()+8*g);"
    "{ OTGEN_HOOK_OUTER(); group<Stream>(in+groupStride*g,out+128*g,plan.field.data()+8*g); }")
  set(S ${SPIN_SOURCE_DIR}/src)
  set_source_files_properties(${CMAKE_BINARY_DIR}/hooked/RouteFast.cpp PROPERTIES INCLUDE_DIRECTORIES "${S}/certified22")
  set_source_files_properties(${CMAKE_BINARY_DIR}/hooked/Fast.cpp PROPERTIES INCLUDE_DIRECTORIES "${S}/rs8")
  # outline=1 keeps the hook work out of SPIN's loops (diagnoses register pressure).
  function(otgen_fused_target name outline)
    add_executable(${name} fused.cpp
      ${S}/Cpu.cpp ${S}/rs8/SetupScalar.cpp ${S}/rs8/Portable.cpp ${S}/rs8/ForwardFast.cpp
      ${S}/certified22/Setup.cpp ${S}/certified22/SetupPortable.cpp ${S}/certified22/Scalar.cpp
      ${S}/certified22/SetupAvx2.cpp ${S}/certified22/SetupFast.cpp
      ${CMAKE_BINARY_DIR}/hooked/RouteFast.cpp ${CMAKE_BINARY_DIR}/hooked/Fast.cpp)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${S} ${SPIN_SOURCE_DIR}/include)
    target_compile_definitions(${name} PRIVATE SPIN_AVX512=1 SPIN_WORKSPACE_ROUTING_OPT=1 OTGEN_HOOK_OUTLINE=${outline})
    target_compile_options(${name} PRIVATE ${OTGEN_ARCH_LIST} -O3 -fno-strict-aliasing)
  endfunction()
  otgen_fused_target(otgen_fused 0)
  otgen_fused_target(otgen_fused_outline 1)
endif()
