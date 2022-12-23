#!/bin/sh
# install plink2 for roundtrip testing
set -ex
cd scripts
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip -P scripts
unzip scripts/plink2_linux_avx2_20221024.zip -d scripts
mv scripts/plink2 /usr/local/bin
# install boost 1.8.0 for tests
wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz -P scripts
tar xf scripts/boost_1_80_0.tar.gz -C scripts
mv scripts/boost_1_80_0 /usr/local/boost
