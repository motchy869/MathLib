cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set(CMAKE_CXX_COMPILER "clang++")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Set the project name and language.
project(MotchyMathLib CXX)

add_subdirectory(include)
add_subdirectory(test)

enable_testing()
add_test(NAME test_analysis COMMAND test_analysis)
add_test(NAME test_linAlg COMMAND test_linAlg)
add_test(NAME test_sigProc COMMAND test_sigProc)
add_test(NAME test_rlsFilter COMMAND test_rlsFilter)
