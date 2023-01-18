cmake_minimum_required(VERSION 3.23)

# This file is used to allow local interactive execution and debugging of the BOOST tests
# for the C++ code from within CLion. The actual build of the JNI and native components for
# this project is managed by the gradle script and dev.nokee plugins.
project(pgen_lib)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/main/headers)
include_directories(src/main/public)
include_directories(/usr/local/boost)

add_executable(pgen_lib
        # headers for the C++ public API (callable by the JNI layer)
        src/main/public/pgenIO.h
        src/main/public/pgenContext.h

        # implementation of the C++ public API (callable by the JNI layer)
        src/main/cpp/pgenIO.cc

        # plink headers
        src/main/headers/pgenlib_ffi_support.h
        src/main/headers/pgenlib_misc.h
        src/main/headers/pgenlib_read.h
        src/main/headers/pgenlib_write.h
        src/main/headers/plink2_base.h
        src/main/headers/plink2_bits.h

        # plink code
        src/main/cpp/pgenlib_ffi_support.cc
        src/main/cpp/pgenlib_misc.cc
        src/main/cpp/pgenlib_read.cc
        src/main/cpp/pgenlib_write.cc
        src/main/cpp/plink2_base.cc
        src/main/cpp/plink2_bits.cc

        # test code
        /usr/local/boost/boost/test/included/unit_test.hpp
        src/test/cpp/test_pgenlib_write.cc src/main/public/pgenException.h)
