cmake_minimum_required(VERSION 3.23)
project(pgen_lib)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/main/headers)
include_directories(/usr/local/boost)

add_executable(pgen_lib
        /usr/local/boost/boost/test/included/unit_test.hpp
        src/main/cpp/pgenlib_ffi_support.cc
        src/main/cpp/pgenlib_misc.cc
        src/main/cpp/pgenlib_read.cc
        src/main/cpp/pgenlib_write.cc
        src/main/cpp/plink2_base.cc
        src/main/cpp/plink2_bits.cc
        src/main/headers/pgenlib_ffi_support.h
        src/main/headers/pgenlib_misc.h
        src/main/headers/pgenlib_read.h
        src/main/headers/pgenlib_write.h
        src/main/headers/plink2_base.h
        src/main/headers/plink2_bits.h
        src/test/cpp/test_pgenlib_write.cc)
