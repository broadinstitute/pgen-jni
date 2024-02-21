# update the local repo with local copies of the plink-ng 
# run this with the root pgen-jni folder as the CWD, assumes plink2-ng clone is parallel with pgen-jni
# headers
cp ../plink-ng/2.0/include/plink2_base.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/plink2_bits.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/pgenlib_misc.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/pgenlib_read.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/include/pgenlib_write.h ./pgen-lib/src/main/headers
cp ../plink-ng/2.0/pgenlib_ffi_support.h ./pgen-lib/src/main/headers

# cpp source
cp ../plink-ng/2.0/include/plink2_base.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/plink2_bits.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/pgenlib_misc.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/pgenlib_read.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/include/pgenlib_write.cc ./pgen-lib/src/main/cpp
cp ../plink-ng/2.0/pgenlib_ffi_support.cc ./pgen-lib/src/main/cpp

# record the commit hash of this update
echo "The plink-ng source in this repo is from plink-ng commit:\n" > README_plink_commit_hash.md
git -C ../plink-ng/ log -1 >> README_plink_commit_hash.md

