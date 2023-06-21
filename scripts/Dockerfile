FROM ubuntu:22.04
LABEL Description="Environment used to build the pgen-jni binary components for Linux"

# NOTE: The versions of plink2 and boost used in here and in install_dependencies.sh
# should be kept in sync.

ENV HOME /root

SHELL ["/bin/bash", "-c"]

# RUN apt update && apt upgrade
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    build-essential \
    clang \
    cmake \
    git \
    openjdk-17-jdk \
    unzip \
    wget
RUN apt-get -y install --reinstall ca-certificates
RUN java --version

# install boost 1.80.0 (required to build) and plink2 executeable (required for tests)
RUN cd ${HOME} && \
    wget --no-check-certificate -q https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz && \
    tar xf boost_1_80_0.tar.gz && \
    mv boost_1_80_0 /usr/local/boost && \
    rm -f boost_1_80_0.tar.gz

# install plink2 executeable (required to run tests)
# use wget -q to supress wget from emitting aws keys into public logs
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha4/plink2_linux_avx2_20230426.zip -P scripts && \
    unzip scripts/plink2_linux_avx2_20230426.zip -d scripts && \
    mv scripts/plink2 /usr/local/bin

# For some reason, the presence of this boost file seems to cause the build of the linux test component (gradle
# task :pgen-lib:compileTestLinuxCpp) to fail. There are numerous online reports of various products failing to
# build on ubuntu due to the presence of the same file. There may be a better solution, but removing it seems to
# work fine for now.
RUN rm -r /usr/local/boost/libs/wave/test/testwave/testfiles/utf8-test-??????????/

