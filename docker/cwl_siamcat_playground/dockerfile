#!/bin/Rscript
FROM ubuntu:16.04

# Install essential unix packages and dependencies:
RUN apt-get update -y --fix-missing
RUN apt-get install -y build-essential \
    wget \
    unzip \
    tar \
    make \
    less \
    nano \
    nodejs \
    npm \
    git
RUN apt-get install -y r-base
# needed for cwltool, toil:
RUN apt-get -y install python-dev \
    python-pip

# Install cwl runners:
ENV TOIL_VERSION 3.13.0
ENV CWLTOOL_VERSION 1.0.20180116213856

RUN pip install toil==${TOIL_VERSION}
RUN pip install cwltool==${CWLTOOL_VERSION}

# Download and install SIAMCAT
ENV SIAMCAT_URL "https://git.embl.de/grp-zeller/SIAMCAT/repository/v0.3.1/archive.tar.gz"
RUN cd /usr/local/ \
    && wget -O siamcat.tar.gz ${SIAMCAT_URL} \
    && tar -xvzf siamcat.tar.gz \
    && mv ./SIAMCAT* ./SIAMCAT \
    && Rscript ./SIAMCAT/Rscript_flavor/00_setup.r ./siamcat.tar.gz \
    && rm siamcat.tar.gz \
    && chmod 700 ./SIAMCAT/Rscript_flavor/*
ENV PATH="/usr/local/SIAMCAT/Rscript_flavor:${PATH}"