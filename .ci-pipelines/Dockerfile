FROM ubuntu:20.04
LABEL maintainer="GEOS-Chem Support Team <geos-chem-support@as.harvard.edu>"
WORKDIR /kpp


RUN apt-get update && apt-get install -y \
    git gcc gfortran make flex bison \
    && rm -rf /var/lib/apt/lists/*
# The last line cleans up cache to reduce container size
# No need to run `apt-get clean` as it is done automatically in official ubuntu image
# See more at https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#run

# To test build locally uncomment COPY . . and run the following commands:
# `docker build -f ./.ci-pipelines/Dockerfile -t kpp-build .` from top level of project
# To exec into the container run: `docker run -it --entrypoint /bin/bash kpp-build`
COPY . .
RUN ls -lha /kpp

# use customized Makefile for ci testing (just switches compiler to gfortran)
RUN rm /kpp/examples/Makefile_saprc_f90 \
    && cp /kpp/.ci-pipelines/pipeline_Makefile_saprc_f90 /kpp/examples/Makefile_saprc_f90

# KPP environment variables
ENV KPP_HOME=/kpp
ENV CC=gcc CXX=g++ FC=gfortran
ENV F90=$FC F9X=$FC F77=$FC

# Build kpp executable update ensure testing script is executable
RUN make && chmod +x /kpp/.ci-pipelines/ci-testing-script.sh 

# Run Tests
RUN /kpp/.ci-pipelines/ci-testing-script.sh

