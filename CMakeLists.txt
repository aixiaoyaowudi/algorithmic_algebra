cmake_minimum_required(VERSION 3.10)

project(Algorithmic_Algebra)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(header)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

find_package(GMP REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY output)


add_executable(test test/test.cpp)
target_link_libraries(test gmp gmpxx)