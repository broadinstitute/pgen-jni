FROM ubuntu:22.04
LABEL Description="Environment used to build the pgen-jni binary components for Linux"

ENV HOME /root

SHELL ["/bin/bash", "-c"]

RUN apt update && apt upgrade
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    build-essential \
    clang \
    cmake \
    git \
    openjdk-17-jdk \
    wget
RUN apt-get -y install --reinstall ca-certificates
RUN java --version

# install boost 1.80.0 for tests
RUN cd ${HOME} && \
    wget --no-check-certificate -q https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz && \
    tar xf boost_1_80_0.tar.gz && \
    mv boost_1_80_0 /usr/local/boost && \
    rm -f boost_1_80_0.tar.gz

# For some reason, the presence of this boost file seems to cause the build of the linux test component (gradle
# task :pgen-lib:compileTestLinuxCpp) to fail. There are numerous online reports of various products failing to
# build on ubuntu due to the presence of the same file. There may be a better solution, but removing it seems to
# work fine for now.
RUN rm -r /usr/local/boost/libs/wave/test/testwave/testfiles/utf8-test-??????????/

