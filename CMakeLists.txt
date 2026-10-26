cmake_minimum_required(VERSION 3.18)
project(fastdiv VERSION 1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 26)

if(NOT MSVC)
  add_compile_options(-Wall -Wextra -O3 -march=native)
endif()

option(SANITIZE "Enable -fsanitize=signed-integer-overflow" OFF)
if(SANITIZE)
  add_compile_options(-fsanitize=signed-integer-overflow)
  add_link_options(-fsanitize=signed-integer-overflow)
endif()

include_directories(${CMAKE_SOURCE_DIR})

# Clang on Windows needs its compiler-rt builtins library for 128-bit integer
# helpers (__udivti3, __umodti3, etc.). The default link line here uses
# -nostdlib, so pull it in explicitly. --print-libgcc-file-name is unreliable
# on Windows (often returns a bare "libgcc.a"), so locate the lib via the
# resource dir and glob for the arch-specific builtins file.
if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -print-resource-dir
    OUTPUT_VARIABLE CLANG_RESOURCE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  file(GLOB CLANG_BUILTINS_LIB "${CLANG_RESOURCE_DIR}/lib/windows/clang_rt.builtins-*.lib")
  if(CLANG_BUILTINS_LIB)
    message(STATUS "Linking clang builtins: ${CLANG_BUILTINS_LIB}")
    link_libraries("${CLANG_BUILTINS_LIB}")
  else()
    message(WARNING "clang_rt.builtins-*.lib not found under ${CLANG_RESOURCE_DIR}/lib/windows")
  endif()
endif()

# The parameter solvers run in constexpr (including 256-bit arithmetic for the
# 64-bit Robison solver), which can exceed the default evaluation budget.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-fconstexpr-steps=500000000)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fconstexpr-ops-limit=2000000000 -fconstexpr-loop-limit=100000000)
elseif(MSVC)
  # /bigobj: div_bench instantiates thousands of benchmark functions.
  add_compile_options(/constexpr:steps500000000 /bigobj /Zc:__cplusplus)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#---------------------------------------------------------------------------
# Google Benchmark
#---------------------------------------------------------------------------

include(FetchContent)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.9.1
  SOURCE_DIR     ${CMAKE_SOURCE_DIR}/tp/benchmark
)
set(BENCHMARK_ENABLE_TESTING OFF)
# Pi Zero and other constrained targets can't run CMake probe executables during
# configure; pre-set the regex backend so benchmark skips its try_run checks.
if(NOT HAVE_POSIX_REGEX AND NOT HAVE_STD_REGEX AND NOT HAVE_GNU_POSIX_REGEX)
  set(HAVE_POSIX_REGEX 1 CACHE INTERNAL "")
endif()
FetchContent_MakeAvailable(benchmark)

#---------------------------------------------------------------------------
# Benchmarks / Tests
#---------------------------------------------------------------------------

add_executable(div_test div_test.cpp)

find_package(Threads REQUIRED)
add_executable(div_report div_report.cpp)
target_link_libraries(div_report Threads::Threads)

add_executable(div_bench div_bench.cpp)
target_link_libraries(div_bench benchmark::benchmark)

# Keep the benchmark honest: never let the toolchain merge identical-looking
# functions or outline repeated sequences, or one algorithm's code could be
# replaced by another's (many divisors produce byte-identical bodies).
if(APPLE)
  target_link_options(div_bench PRIVATE -Wl,-no_deduplicate)
elseif(MSVC)
  target_link_options(div_bench PRIVATE /OPT:NOICF)
else()
  # --icf is a gold/lld extension; GNU ld (bfd) rejects it outright.
  include(CheckLinkerFlag)
  check_linker_flag(CXX "-Wl,--icf=none" LINKER_SUPPORTS_ICF_NONE)
  if(LINKER_SUPPORTS_ICF_NONE)
    target_link_options(div_bench PRIVATE -Wl,--icf=none)
  endif()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    target_compile_options(div_bench PRIVATE -mno-outline)
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(div_bench PRIVATE -fno-ipa-icf)
endif()
