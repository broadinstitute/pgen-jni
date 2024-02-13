#!/bin/sh

# NOTE: The versions of plink2 and boost used here should be kept in sync with the parallel srei for linux (install_linux_dependencies.sh).

# install plink2 for roundtrip testing
set -ex
cd scripts
# use wget -q to supress wget from emitting aws keys into public logs
wget -q https://s3.amazonaws.com/plink2-assets/plink2_mac_20240205.zip -P scripts
unzip scripts/plink2_mac_20240205.zip -d scripts
sudo mv scripts/plink2 /usr/local/bin

# install boost 1.8.0 for tests
wget -q https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz -P scripts
tar xf scripts/boost_1_80_0.tar.gz -C scripts
sudo mv scripts/boost_1_80_0 /usr/local/boost