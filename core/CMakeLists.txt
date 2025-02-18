cmake_minimum_required(VERSION 3.26)
project(core)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



#include_directories("D:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Include")
#link_directories("D:\\Program Files (x86)\\Microsoft SDKs\\MPI\\Lib\\x64")
#link_libraries(msmpi)s
include_directories(${PROJECT_BINARY_DIR}/../include)
file(GLOB src_list "${PROJECT_BINARY_DIR}/../src/*.cpp" "${PROJECT_BINARY_DIR}/../src/*.cc")
add_executable(main ${src_list})
target_link_libraries(main OpenMP::OpenMP_CXX)