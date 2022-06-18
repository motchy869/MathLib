cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Include and link files directory (for googletest)
# You may have to change path according to your environment.
if (UNIX)
    set(CMAKE_GOOGLE_TEST_INCLUDE_DIR "/usr/local/include/")
    set(CMAKE_GTEST_LIB_DIR "/usr/local/lib/gtest/")
    set(CMAKE_GMOCK_LIB_DIR "/usr/local/lib/gmock/")
elseif (WIN32)
    set(CMAKE_GOOGLE_TEST_INCLUDE_DIR "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/include/")
    set(CMAKE_GTEST_LIB_DIR "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gtest/")
    set(CMAKE_GMOCK_LIB_DIR "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gmock/")
endif ()
include_directories(${CMAKE_GOOGLE_TEST_INCLUDE_DIR})
link_directories(${CMAKE_GTEST_LIB_DIR})
link_directories(${CMAKE_GMOCK_LIB_DIR})

# Generate executable file
add_executable(test_analysis ${CMAKE_CURRENT_SOURCE_DIR}/test_analysis.cpp)
add_executable(test_linAlg ${CMAKE_CURRENT_SOURCE_DIR}/test_linAlg.cpp)
add_executable(test_sigProc ${CMAKE_CURRENT_SOURCE_DIR}/test_sigProc.cpp)
add_executable(test_RlsFilter ${CMAKE_CURRENT_SOURCE_DIR}/test_RlsFilter.cpp)

# Link libraries
target_link_libraries(test_analysis motchyMathLib gtest gtest_main pthread)
target_link_libraries(test_linAlg motchyMathLib gtest gtest_main pthread)
target_link_libraries(test_sigProc motchyMathLib gtest gtest_main pthread)
target_link_libraries(test_RlsFilter motchyMathLib gtest gtest_main pthread)
