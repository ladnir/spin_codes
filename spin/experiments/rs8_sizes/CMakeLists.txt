cmake_minimum_required(VERSION 3.20)
project(spin_rs8_sizes LANGUAGES CXX)
find_package(spin CONFIG REQUIRED)
add_executable(spin_rs8_sizes bench.cpp Setup15.cpp Route15.cpp
  ../k20_two_ms/rs8/SetupScalar.cpp ../k20_two_ms/rs8/Fast.cpp)
target_link_libraries(spin_rs8_sizes PRIVATE spin::spin)
target_compile_features(spin_rs8_sizes PRIVATE cxx_std_20)
target_compile_options(spin_rs8_sizes PRIVATE -mavx2 -mavx512f -mavx512vl
  -mavx512bw -mavx512dq -mavx512vbmi -mgfni -mtune=znver4
  -fno-strict-aliasing -fno-lto)
enable_testing()
foreach(mode IN ITEMS baseline15 rs8-15 baseline22 rs8-22)
  foreach(size IN ITEMS 2176 6144 65536 262144)
    foreach(seed IN ITEMS 17 41)
      add_test(NAME ${mode}_${size}_${seed}
        COMMAND spin_rs8_sizes ${mode} ${size} ${seed} 0)
      set_tests_properties(${mode}_${size}_${seed} PROPERTIES RUN_SERIAL TRUE SKIP_RETURN_CODE 77)
    endforeach()
  endforeach()
endforeach()
