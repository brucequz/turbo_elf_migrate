cmake_minimum_required(VERSION 3.19.1)
project(turbo-elf VERSION 0.1.0 LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 14)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -fsanitize=address")

file(GLOB SOURCES "src/*.cpp")
include_directories(${CMAKE_SOURCE_DIR}/include/*)

add_executable(turbo ${SOURCES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
