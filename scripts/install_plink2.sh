#!/bin/sh
# install plink2 for roundtrip testing
set -ex
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip
unzip plink2_linux_avx2_20221024.zip
mv plink2 /usr/local/bin
