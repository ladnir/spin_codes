cmake_minimum_required(VERSION 3.20)

project(GolayBa3Rm2SubB240K20 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SPIN_ARCH "znver4" CACHE STRING "GCC/Clang -march value")
set(STRUCTURED_SPIN_IMPLEMENTATION_ROOT
    "${CMAKE_CURRENT_LIST_DIR}/../../../constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation"
    CACHE PATH "Path to the cleaned Structured SPIN implementation")

add_library(golay_ba3_rm2sub_b240_k20 STATIC
    libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.cpp)
target_include_directories(golay_ba3_rm2sub_b240_k20 PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${STRUCTURED_SPIN_IMPLEMENTATION_ROOT}"
    "${STRUCTURED_SPIN_IMPLEMENTATION_ROOT}/vendor")
target_compile_features(golay_ba3_rm2sub_b240_k20 PUBLIC cxx_std_20)

add_library(structured_spin_b256_baseline STATIC
    "${STRUCTURED_SPIN_IMPLEMENTATION_ROOT}/libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.cpp")
target_include_directories(structured_spin_b256_baseline PUBLIC
    "${STRUCTURED_SPIN_IMPLEMENTATION_ROOT}"
    "${STRUCTURED_SPIN_IMPLEMENTATION_ROOT}/vendor")
target_compile_features(structured_spin_b256_baseline PUBLIC cxx_std_20)

add_executable(golay_ba3_correctness tests/golay_ba3_correctness.cpp)
target_link_libraries(golay_ba3_correctness PRIVATE golay_ba3_rm2sub_b240_k20)

add_executable(golay_ba3_public_header_smoke tests/public_header_smoke.cpp)
target_link_libraries(golay_ba3_public_header_smoke PRIVATE golay_ba3_rm2sub_b240_k20)

add_executable(golay_ba3_comparison benchmarks/golay_ba3_comparison.cpp)
target_link_libraries(golay_ba3_comparison PRIVATE
    golay_ba3_rm2sub_b240_k20 structured_spin_b256_baseline)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(target golay_ba3_rm2sub_b240_k20 structured_spin_b256_baseline
            golay_ba3_correctness golay_ba3_comparison)
        target_compile_options(${target} PRIVATE
            "-march=${SPIN_ARCH}" -mavx2 -mvpclmulqdq -mpclmul)
    endforeach()
elseif(MSVC)
    foreach(target golay_ba3_rm2sub_b240_k20 structured_spin_b256_baseline
            golay_ba3_correctness golay_ba3_comparison)
        target_compile_options(${target} PRIVATE /arch:AVX2)
    endforeach()
endif()

if(UNIX)
    target_link_libraries(golay_ba3_correctness PRIVATE pthread)
    target_link_libraries(golay_ba3_comparison PRIVATE pthread)
endif()

enable_testing()
add_test(NAME golay_ba3_correctness COMMAND golay_ba3_correctness)
add_test(NAME golay_ba3_public_header_smoke COMMAND golay_ba3_public_header_smoke)
