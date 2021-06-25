cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Include files directory (for gtest)
include_directories(/usr/local/include/gtest/)
include_directories(/usr/local/include/gmock/)

# Library directory (for gtest)
link_directories(/usr/local/lib/gtest/)
link_directories(/usr/local/lib/gmock/)

# Generate executable file
add_executable(test_linAlg ${CMAKE_CURRENT_SOURCE_DIR}/test_linAlg.cpp)
add_executable(test_sigProc ${CMAKE_CURRENT_SOURCE_DIR}/test_sigProc.cpp)
add_executable(test_rlsFilter ${CMAKE_CURRENT_SOURCE_DIR}/test_rlsFilter.cpp)

# Link libraries
target_link_libraries(test_linAlg gtest gtest_main pthread)
target_link_libraries(test_sigProc gtest gtest_main pthread)
target_link_libraries(test_rlsFilter gtest gtest_main pthread)
