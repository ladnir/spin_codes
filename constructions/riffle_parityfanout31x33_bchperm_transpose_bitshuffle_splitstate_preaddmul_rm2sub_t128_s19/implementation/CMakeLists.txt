cmake_minimum_required(VERSION 3.20)

project(StructuredSpinB256T128S19 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SPIN_ARCH "znver4" CACHE STRING "GCC/Clang -march value for Structured SPIN")

add_library(structured_spin_b256_t128_s19 STATIC
    libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.cpp)
target_include_directories(structured_spin_b256_t128_s19 PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vendor")
target_compile_features(structured_spin_b256_t128_s19 PUBLIC cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(structured_spin_b256_t128_s19 PRIVATE
        "-march=${SPIN_ARCH}"
        -mavx2
        -mvpclmulqdq
        -mpclmul)
elseif(MSVC)
    target_compile_options(structured_spin_b256_t128_s19 PRIVATE /arch:AVX2)
endif()

add_executable(structured_spin_correctness
    tests/structured_spin_correctness.cpp)
target_link_libraries(structured_spin_correctness PRIVATE
    structured_spin_b256_t128_s19)

add_executable(structured_spin_benchmark
    benchmarks/structured_spin_benchmark.cpp)
target_link_libraries(structured_spin_benchmark PRIVATE
    structured_spin_b256_t128_s19)

add_executable(structured_spin_public_header_smoke
    tests/public_header_smoke.cpp)
target_link_libraries(structured_spin_public_header_smoke PRIVATE
    structured_spin_b256_t128_s19)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(target structured_spin_correctness structured_spin_benchmark)
        target_compile_options(${target} PRIVATE
            "-march=${SPIN_ARCH}"
            -mavx2
            -mvpclmulqdq
            -mpclmul)
    endforeach()
elseif(MSVC)
    target_compile_options(structured_spin_correctness PRIVATE /arch:AVX2)
    target_compile_options(structured_spin_benchmark PRIVATE /arch:AVX2)
endif()

if(UNIX)
    target_link_libraries(structured_spin_correctness PRIVATE pthread)
    target_link_libraries(structured_spin_benchmark PRIVATE pthread)
endif()

enable_testing()
add_test(NAME structured_spin_correctness COMMAND structured_spin_correctness)
add_test(NAME structured_spin_public_header_smoke
    COMMAND structured_spin_public_header_smoke)
