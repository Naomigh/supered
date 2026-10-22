cmake_minimum_required(VERSION 3.20)
project(supered LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SUPERED_LEAP_SOURCES
    leap_origin/SIMD_ED.cc
    leap_origin/SHD.cc
    leap_origin/bit_convert.c
    leap_origin/shift.c
    leap_origin/mask.c
    leap_origin/popcount.c
    leap_origin/print.c)
# LEAP's historical .c files contain C++ intrinsics.
set_source_files_properties(${SUPERED_LEAP_SOURCES} PROPERTIES LANGUAGE CXX)

add_library(supered STATIC ${SUPERED_LEAP_SOURCES})
add_library(supered::supered ALIAS supered)
target_compile_features(supered PUBLIC cxx_std_20)
target_include_directories(supered PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/leap_origin)
target_compile_options(supered PUBLIC -march=native PRIVATE
    -O3 -Wall -Wextra -Wno-unused-function -Wno-unused-variable -Wno-sign-compare)
set_property(TARGET supered PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    option(SUPERED_BUILD_EXAMPLES "Build the basic usage example" ON)
else()
    option(SUPERED_BUILD_EXAMPLES "Build the basic usage example" OFF)
endif()

if(SUPERED_BUILD_EXAMPLES)
    add_executable(supered_example examples/basic.cpp)
    target_link_libraries(supered_example PRIVATE supered::supered)
    target_compile_options(supered_example PRIVATE -O3)
    set_property(TARGET supered_example PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
