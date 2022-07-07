# This dockerfile is used to run the CI unit tests on azure pipelines.
# To test build locally run the following commands:
# `docker build -f ./.ci-pipelines/Dockerfile --platform=amd64 -t kpp-build .` from top level of project
# To exec into the container run: `docker run -it --entrypoint /bin/bash kpp-build`

# Set OS, a reference to maintainers, and the working directory.
#
# NOTE: Add --platform=amd64 (in .build-testing.yml and in your local
# aliases for 'docker build'.  This will make sure that the C-I tests always
# run on AMD64 hardware.  Otherwise, Docker will use the architecture of
# your current hardware, which might cause differences in library paths
# (e.g. /lib/aarch64-linux-gnu on Mac/M1 vs. /lib/x86_64-linux-gnu on Linux).
#   -- Lucas Estrada & Bob Yantosca (06 Jul 2022)
FROM ubuntu:20.04
LABEL maintainer="GEOS-Chem Support Team <geos-chem-support@g.harvard.edu>"
WORKDIR /kpp

# Install KPP software dependencies
RUN apt-get update && apt-get install -y git gcc gfortran make flex bison

COPY . .

# Set environment variables for building KPP
ENV KPP_HOME=/kpp
ENV CC=gcc CXX=g++ FC=gfortran
ENV F90=$FC F9X=$FC F77=$FC
ENV KPP_FLEX_LIB_DIR=/lib/x86_64-linux-gnu

# For C-I tests MAX_EQN and MAX_SPECIES must be < 1024:
# as this will fit within the memory limits on Azure.
RUN sed -i 's/#define MAX_EQN .*/#define MAX_EQN 1023/g' /kpp/src/gdata.h \
    && sed -i 's/#define MAX_SPECIES .*/#define MAX_SPECIES 1023/g' /kpp/src/gdata.h

# Build KPP executable and ensure testing scripts are executable
RUN cd /kpp/src/ && make
RUN chmod +x /kpp/.ci-pipelines/ci-common-defs.sh
RUN chmod +x /kpp/.ci-pipelines/ci-testing-script.sh

# Run C-I tests
RUN /kpp/.ci-pipelines/ci-testing-script.sh

