# update the local repo with local copies of the plink-ng 
# run this with the root pgen-jni folder as the CWD, assumes plink2-ng clone is parallel with pgen-jni
# headers
cp ../plink-ng/2.0/include/plink2_base.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/plink2_bits.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/pgenlib_misc.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/pgenlib_write.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/pgenlib_ffi_support.h ./pgen-lib/src/main/headers

# cpp source
cp ../plink-ng/2.0/include/plink2_base.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/plink2_bits.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/pgenlib_misc.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/pgenlib_read.cc ./pgen-lib/src/main/cpp/pgenlib_read.cc
cp ../plink-ng/2.0/include/pgenlib_write.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/pgenlib_ffi_support.cc ./pgen-lib/src/main/cpp/pgenlib_ffi_support.cc

