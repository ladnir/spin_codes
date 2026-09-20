# Initial cache for a fresh optional forward/wide build. No global defaults change.
# Supply SPIN_BIDIRECTIONAL_SOURCE on the command line; upstream is copied, not edited.
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set(SPIN_BCH_AVX512 ON CACHE BOOL "Compile runtime-dispatched four-row BCH")
set(SPIN_BUILD_WIDE ON CACHE BOOL "Build 256/512-bit forward kernels")
set(SPIN_FORWARD_FOUR ON CACHE BOOL "Four-row 128-bit forward BCH")
set(SPIN_FORWARD_DIRECT ON CACHE BOOL "Direct forward routing through K=2^18")
set(SPIN_FORWARD_SCHEDULE dfs CACHE STRING "128-bit forward BCH schedule")
set(SPIN_WIDE_SCHEDULE dfs CACHE STRING "Wide forward BCH schedule")
