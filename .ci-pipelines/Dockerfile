# This dockerfile is used to run the CI unit tests on azure pipelines.
# To test build locally run the following commands:
# `docker build -f ./.ci-pipelines/Dockerfile -t kpp-build .` from top level of project
# To exec into the container run: `docker run -it --entrypoint /bin/bash kpp-build`
FROM ubuntu:20.04
LABEL maintainer="GEOS-Chem Support Team <geos-chem-support@as.harvard.edu>"
WORKDIR /kpp

# install software pre reqs
RUN apt-get update && apt-get install -y \
    git gcc gfortran make flex bison

COPY . .

# KPP environment variables
ENV KPP_HOME=/kpp
ENV CC=gcc CXX=g++ FC=gfortran
ENV F90=$FC F9X=$FC F77=$FC

# for CI, MAX_EQN and MAX_SPECIES must be < 1024:
RUN sed -i 's/#define MAX_EQN .*/#define MAX_EQN 1023/g' /kpp/src/gdata.h \
    && sed -i 's/#define MAX_SPECIES .*/#define MAX_SPECIES 1023/g' /kpp/src/gdata.h

# switch compiler to use gfortran instead of ifort
RUN sed -i 's/COMPILER = INTEL/#COMPILER = INTEL/g' /kpp/util/Makefile.f90 \
    && sed -i 's/#COMPILER = GFORTRAN/COMPILER = GFORTRAN/g' /kpp/util/Makefile.f90

# Build kpp executable update ensure testing script is executable
RUN make && chmod +x /kpp/.ci-pipelines/ci-testing-script.sh 

# Run Tests
RUN /kpp/.ci-pipelines/ci-testing-script.sh

