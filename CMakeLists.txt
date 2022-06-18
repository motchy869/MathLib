cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if (UNIX)
    set(CMAKE_CXX_COMPILER "clang++")
    string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address")
elseif (WIN32)
    set(CMAKE_CXX_COMPILER "g++")
else ()
    message(FATAL_ERROR "Unknown platform.")
endif ()

set(CMAKE_CXX_STANDARD 17)
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Wfatal-errors -Winline -Wold-style-cast -Wredundant-decls -Wshadow -Wswitch-default -Wswitch-enum -Wundef -Wunused-variable")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g -fstack-protector -ftrapv")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O3 -DNDEBUG")

# Set the project name and language.
project(MathLib CXX)

add_subdirectory(include)
add_subdirectory(test)

enable_testing()
add_test(NAME test_analysis COMMAND test_analysis)
add_test(NAME test_linAlg COMMAND test_linAlg)
add_test(NAME test_sigProc COMMAND test_sigProc)
add_test(NAME test_RlsFilter COMMAND test_RlsFilter)
