cmake_minimum_required(VERSION 3.20)
project(spin_char2_field LANGUAGES CXX)
option(SPIN_CHAR2_CLMUL_REDUCTION "Use PCLMUL/VPCLMUL for the native GF128 overflow fold" OFF)
option(SPIN_CHAR2_SCALAR_LOOKUP "Use scalar affine-table loads instead of AVX512 gathers" OFF)

add_library(spin_char2_field STATIC Portable.cpp Kernel128.cpp Kernel512.cpp)
target_compile_features(spin_char2_field PUBLIC cxx_std_20)
target_include_directories(spin_char2_field PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(spin_char2_field PUBLIC
  SPIN_CHAR2_CLMUL_REDUCTION=$<BOOL:${SPIN_CHAR2_CLMUL_REDUCTION}>
  SPIN_CHAR2_SCALAR_LOOKUP=$<BOOL:${SPIN_CHAR2_SCALAR_LOOKUP}>)
if(MSVC)
  # MSVC lacks per-feature GFNI flags. Kernel128's public availability predicate
  # requires AVX2/GFNI; Kernel512 separately requires OS-enabled AVX512 state.
  set_source_files_properties(Kernel128.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
  set_source_files_properties(Kernel512.cpp PROPERTIES COMPILE_OPTIONS "/arch:AVX512")
  target_compile_options(spin_char2_field PRIVATE /O2 /W4)
else()
  set_source_files_properties(Kernel128.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mgfni")
  set_source_files_properties(Kernel512.cpp PROPERTIES COMPILE_OPTIONS "-mavx2;-mavx512f;-mavx512bw;-mavx512vl;-mgfni")
  target_compile_options(spin_char2_field PRIVATE -O3 -Wall -Wextra -Wpedantic -fno-lto)
  if(SPIN_CHAR2_CLMUL_REDUCTION)
    set_property(SOURCE Kernel128.cpp APPEND PROPERTY COMPILE_OPTIONS -mpclmul)
    set_property(SOURCE Kernel512.cpp APPEND PROPERTY COMPILE_OPTIONS -mpclmul -mvpclmulqdq)
  endif()
endif()

add_executable(spin_char2_field_frontend frontend.cpp)
target_link_libraries(spin_char2_field_frontend PRIVATE spin_char2_field)
enable_testing()
add_test(NAME char2_small_coefficients COMMAND spin_char2_field_frontend --check)
set_tests_properties(char2_small_coefficients PROPERTIES RUN_SERIAL TRUE)

add_executable(spin_char2_full full.cpp)
target_link_libraries(spin_char2_full PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_full PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_full PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
  if(SPIN_CHAR2_CLMUL_REDUCTION)
    target_compile_options(spin_char2_full PRIVATE -mpclmul -mvpclmulqdq)
  endif()
endif()
add_test(NAME char2_full_encoder COMMAND spin_char2_full --check)
set_tests_properties(char2_full_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_batched batched.cpp)
target_link_libraries(spin_char2_batched PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_batched PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_batched PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
  if(SPIN_CHAR2_CLMUL_REDUCTION)
    target_compile_options(spin_char2_batched PRIVATE -mpclmul -mvpclmulqdq)
  endif()
endif()
add_test(NAME char2_batched_encoder COMMAND spin_char2_batched --check)
set_tests_properties(char2_batched_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_subfield subfield.cpp)
target_link_libraries(spin_char2_subfield PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_subfield PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_subfield PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_subfield_encoder COMMAND spin_char2_subfield --check)
set_tests_properties(char2_subfield_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_one_route oneroute.cpp)
target_link_libraries(spin_char2_one_route PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_one_route PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_one_route PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_one_route_encoder COMMAND spin_char2_one_route --check)
set_tests_properties(char2_one_route_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_subfield16 subfield16.cpp)
target_link_libraries(spin_char2_subfield16 PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_subfield16 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_subfield16 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_subfield16_encoder COMMAND spin_char2_subfield16 --check)
set_tests_properties(char2_subfield16_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Separate research candidate: the same GF16 route with a six-symbol state.
add_executable(spin_char2_subfield_wide wide.cpp)
target_link_libraries(spin_char2_subfield_wide PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_subfield_wide PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_subfield_wide PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_subfield_wide_encoder COMMAND spin_char2_subfield_wide --check)
set_tests_properties(char2_subfield_wide_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Independently sampled Toeplitz[64,32] outer; one regional route and GL3 inner.
add_executable(spin_char2_toeplitz toeplitz.cpp)
target_link_libraries(spin_char2_toeplitz PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_toeplitz PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_toeplitz PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_toeplitz_encoder COMMAND spin_char2_toeplitz --check)
set_tests_properties(char2_toeplitz_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Higher-rate research checkpoint; production code is only an optional control.
set(SPIN_CHAR2_PUBLIC_ROOT "" CACHE PATH "Optional standalone SPIN source for matched public-API comparison")
add_executable(spin_char2_cauchy256 cauchy256.cpp)
target_link_libraries(spin_char2_cauchy256 PRIVATE spin_char2_field)
if(SPIN_CHAR2_PUBLIC_ROOT)
  add_subdirectory(${SPIN_CHAR2_PUBLIC_ROOT} ${CMAKE_CURRENT_BINARY_DIR}/public-spin EXCLUDE_FROM_ALL)
  target_link_libraries(spin_char2_cauchy256 PRIVATE spin::spin)
  target_compile_definitions(spin_char2_cauchy256 PRIVATE SPIN_CAUCHY_PUBLIC=1)
endif()
if(MSVC)
  target_compile_options(spin_char2_cauchy256 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_cauchy256 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_cauchy256_encoder COMMAND spin_char2_cauchy256 --check)
set_tests_properties(char2_cauchy256_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Smaller higher-rate outer with an explicit four-row-interleaved source order.
add_executable(spin_char2_cauchy128 cauchy128.cpp)
target_link_libraries(spin_char2_cauchy128 PRIVATE spin_char2_field)
if(SPIN_CHAR2_PUBLIC_ROOT)
  target_link_libraries(spin_char2_cauchy128 PRIVATE spin::spin)
  target_compile_definitions(spin_char2_cauchy128 PRIVATE SPIN_CAUCHY_PUBLIC=1)
endif()
if(MSVC)
  target_compile_options(spin_char2_cauchy128 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_cauchy128 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_cauchy128_encoder COMMAND spin_char2_cauchy128 --check)
set_tests_properties(char2_cauchy128_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Isolated medium-field arithmetic controls; not a full encoder benchmark.
add_executable(spin_char2_medium_field medium_field.cpp)
target_compile_features(spin_char2_medium_field PRIVATE cxx_std_17)
if(MSVC)
  target_compile_options(spin_char2_medium_field PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_medium_field PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mvpclmulqdq -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_medium_field_arithmetic COMMAND spin_char2_medium_field --check)
set_tests_properties(char2_medium_field_arithmetic PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Same certified outer image in an additive-polynomial source basis.
add_executable(spin_char2_cauchy128_coefficients cauchy128_coefficients.cpp)
target_link_libraries(spin_char2_cauchy128_coefficients PRIVATE spin_char2_field)
if(SPIN_CHAR2_PUBLIC_ROOT)
  target_link_libraries(spin_char2_cauchy128_coefficients PRIVATE spin::spin)
  target_compile_definitions(spin_char2_cauchy128_coefficients PRIVATE SPIN_CAUCHY_PUBLIC=1)
endif()
if(MSVC)
  target_compile_options(spin_char2_cauchy128_coefficients PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_cauchy128_coefficients PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_cauchy128_coefficient_encoder COMMAND spin_char2_cauchy128_coefficients --check)
set_tests_properties(char2_cauchy128_coefficient_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Certified K20 research candidate: GRS[128,83] and repaired cubic width8/state4.
add_executable(spin_char2_cubic8 cubic8.cpp)
target_link_libraries(spin_char2_cubic8 PRIVATE spin_char2_field)
if(SPIN_CHAR2_PUBLIC_ROOT)
  target_link_libraries(spin_char2_cubic8 PRIVATE spin::spin)
  target_compile_definitions(spin_char2_cubic8 PRIVATE SPIN_CAUCHY_PUBLIC=1)
endif()
if(MSVC)
  target_compile_options(spin_char2_cubic8 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_cubic8 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_cubic8_encoder COMMAND spin_char2_cubic8 --check)
set_tests_properties(char2_cubic8_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Isolated nominal-rate3/4 research candidate; no production code changes.
add_executable(spin_char2_punctured8 punctured8.cpp)
target_link_libraries(spin_char2_punctured8 PRIVATE spin_char2_field)
option(SPIN_P8_FORCE_INLINE "Inline the four-packet hook to eliminate its scratch roundtrip" ON)
target_compile_definitions(spin_char2_punctured8 PRIVATE SPIN_P8_FORCE_INLINE=$<BOOL:${SPIN_P8_FORCE_INLINE}>)
if(MSVC)
  target_compile_options(spin_char2_punctured8 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_punctured8 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_punctured8_encoder COMMAND spin_char2_punctured8 --check)
set_tests_properties(char2_punctured8_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Separately certified high-rate punctures; frozen 75% target stays unchanged.
add_executable(spin_char2_rate_climb8 rate_climb8.cpp)
target_link_libraries(spin_char2_rate_climb8 PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_rate_climb8 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_rate_climb8 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_rate_climb8_encoder COMMAND spin_char2_rate_climb8 --check)
set_tests_properties(char2_rate_climb8_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_rate_speed8 rate_speed8.cpp)
target_link_libraries(spin_char2_rate_speed8 PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_rate_speed8 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_rate_speed8 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_rate_speed8_encoder COMMAND spin_char2_rate_speed8 --check)
set_tests_properties(char2_rate_speed8_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_rate_speed4 rate_speed4.cpp)
target_link_libraries(spin_char2_rate_speed4 PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_rate_speed4 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_rate_speed4 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_rate_speed4_encoder COMMAND spin_char2_rate_speed4 --check)
set_tests_properties(char2_rate_speed4_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_rate_speed4_next rate_speed4_next.cpp)
target_link_libraries(spin_char2_rate_speed4_next PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_rate_speed4_next PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_rate_speed4_next PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_rate_speed4_next_encoder COMMAND spin_char2_rate_speed4_next --check)
set_tests_properties(char2_rate_speed4_next_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_rate_speed4_profile rate_speed4_profile.cpp)
target_link_libraries(spin_char2_rate_speed4_profile PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_rate_speed4_profile PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_rate_speed4_profile PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_rate_speed4_profile_encoder COMMAND spin_char2_rate_speed4_profile --check)
set_tests_properties(char2_rate_speed4_profile_encoder PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Isolated arithmetic control for the high-distance field exploration.
add_executable(spin_char2_wide_diagonal16 wide_diagonal16.cpp)
target_link_libraries(spin_char2_wide_diagonal16 PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_wide_diagonal16 PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_wide_diagonal16 PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_wide_diagonal16 COMMAND spin_char2_wide_diagonal16 --check)
set_tests_properties(char2_wide_diagonal16 PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Full rate-one-half wide-field cost prototype, not a distance certificate.
add_executable(spin_char2_high_distance_wide high_distance_wide.cpp)
target_link_libraries(spin_char2_high_distance_wide PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_high_distance_wide PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_high_distance_wide PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_high_distance_wide COMMAND spin_char2_high_distance_wide --check)
set_tests_properties(char2_high_distance_wide PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Direction-matched forward control for the high-distance research profile.
add_executable(spin_char2_hd_fast_forward hd_fast_forward.cpp)
target_link_libraries(spin_char2_hd_fast_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_fast_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_fast_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_fast_forward COMMAND spin_char2_hd_fast_forward --check)
set_tests_properties(char2_hd_fast_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_hd_striped_baa hd_striped_baa.cpp)
target_link_libraries(spin_char2_hd_striped_baa PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_striped_baa PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_striped_baa PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_striped_baa COMMAND spin_char2_hd_striped_baa --check)
set_tests_properties(char2_hd_striped_baa PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_hd_affine_two_state hd_affine_two_state.cpp)
target_link_libraries(spin_char2_hd_affine_two_state PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_affine_two_state PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_affine_two_state PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_affine_two_state COMMAND spin_char2_hd_affine_two_state --check)
set_tests_properties(char2_hd_affine_two_state PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Direct SPIN restart: one route and a batched stateful inner, no pre-prefix.
add_executable(spin_char2_hd_blockwise_forward hd_blockwise_forward.cpp)
target_link_libraries(spin_char2_hd_blockwise_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_blockwise_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_blockwise_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_blockwise_forward COMMAND spin_char2_hd_blockwise_forward --check)
set_tests_properties(char2_hd_blockwise_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_hd_blockwise_small_forward hd_blockwise_small_forward.cpp)
target_link_libraries(spin_char2_hd_blockwise_small_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_blockwise_small_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_blockwise_small_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_blockwise_small_forward COMMAND spin_char2_hd_blockwise_small_forward --check)
set_tests_properties(char2_hd_blockwise_small_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Exact beta64 shared cube/FFT map controls; no distance certificate.
add_executable(spin_char2_hd_blockwise_shared_forward hd_blockwise_shared_forward.cpp)
target_link_libraries(spin_char2_hd_blockwise_shared_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_blockwise_shared_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_blockwise_shared_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_blockwise_shared_forward COMMAND spin_char2_hd_blockwise_shared_forward --check)
set_tests_properties(char2_hd_blockwise_shared_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Independent nonzero GF65536 outer-output scalars, fused into region staging.
# Isolated prototype: the existing shared driver and defaults are preserved.
add_executable(spin_char2_hd_blockwise_randomized_forward hd_blockwise_randomized_forward.cpp)
target_link_libraries(spin_char2_hd_blockwise_randomized_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_blockwise_randomized_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_blockwise_randomized_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_blockwise_randomized_forward COMMAND spin_char2_hd_blockwise_randomized_forward --check)
set_tests_properties(char2_hd_blockwise_randomized_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

add_executable(spin_char2_hd_sumrank_forward hd_sumrank_forward.cpp)
target_link_libraries(spin_char2_hd_sumrank_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_sumrank_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_sumrank_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_sumrank_forward COMMAND spin_char2_hd_sumrank_forward --check)
set_tests_properties(char2_hd_sumrank_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Wider systematic-MDS research candidates; separate from the frozen s16 path.
add_executable(spin_char2_hd_wide_mds_forward hd_wide_mds_forward.cpp)
target_link_libraries(spin_char2_hd_wide_mds_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_wide_mds_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_wide_mds_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_wide_mds_forward COMMAND spin_char2_hd_wide_mds_forward --check)
set_tests_properties(char2_hd_wide_mds_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# 35%-distance co-design candidate, isolated from the frozen 40% controls.
add_executable(spin_char2_hd35_forward hd35_forward.cpp)
target_link_libraries(spin_char2_hd35_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd35_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd35_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd35_forward COMMAND spin_char2_hd35_forward --check)
set_tests_properties(char2_hd35_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)

# Grouped small RS rows with independent cross-row column randomizers.
add_executable(spin_char2_hd_grouped_forward hd_grouped_forward.cpp)
target_link_libraries(spin_char2_hd_grouped_forward PRIVATE spin_char2_field)
if(MSVC)
  target_compile_options(spin_char2_hd_grouped_forward PRIVATE /arch:AVX512 /O2)
else()
  target_compile_options(spin_char2_hd_grouped_forward PRIVATE -O3 -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto)
endif()
add_test(NAME char2_hd_grouped_forward COMMAND spin_char2_hd_grouped_forward --check)
set_tests_properties(char2_hd_grouped_forward PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
