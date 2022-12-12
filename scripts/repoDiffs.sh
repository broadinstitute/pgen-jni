# diff the local copies of the plink-ng source with the real versions in the repo
# for use when refreshing the sources
# headers
diff ./pgen-lib/src/main/headers/plink2_base.h ../plink-ng/2.0/include/plink2_base.h
diff ./pgen-lib/src/main/headers/plink2_bits.h ../plink-ng/2.0/include/plink2_bits.h
diff ./pgen-lib/src/main/headers/pgenlib_misc.h ../plink-ng/2.0/include/pgenlib_misc.h
diff ./pgen-lib/src/main/headers/pgenlib_read.h ../plink-ng/2.0/include/pgenlib_read.h
diff ./pgen-lib/src/main/headers/pgenlib_write.h ../plink-ng/2.0/include/pgenlib_write.h
diff ./pgen-lib/src/main/headers/pgenlib_ffi_support.h ../plink-ng/2.0/pgenlib_ffi_support.h

# cpp source
diff ./pgen-lib/src/main/cpp/plink2_base.cc ../plink-ng/2.0/include/plink2_base.cc
diff ./pgen-lib/src/main/cpp/plink2_bits.cc ../plink-ng/2.0/include/plink2_bits.cc
diff ./pgen-lib/src/main/cpp/pgenlib_misc.cc ../plink-ng/2.0/include/pgenlib_misc.cc
diff ./pgen-lib/src/main/cpp/pgenlib_read.cc ../plink-ng/2.0/include/pgenlib_read.cc
diff ./pgen-lib/src/main/cpp/pgenlib_write.cc ../plink-ng/2.0/include/pgenlib_write.cc
diff ./pgen-lib/src/main/cpp/pgenlib_ffi_support.cc ../plink-ng/2.0/pgenlib_ffi_support.cc

