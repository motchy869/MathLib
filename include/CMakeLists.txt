cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

add_library(motchyMathLib INTERFACE)
target_include_directories(motchyMathLib INTERFACE ${PROJECT_SOURCE_DIR}/include)
