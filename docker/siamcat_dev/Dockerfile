#!/bin/Rscript
FROM ubuntu:16.04

# Install essential unix packages and dependencies:
RUN apt-get update -y --fix-missing
RUN apt-get install -y build-essential \
    wget \
    tar
RUN apt-get install -y r-base


# Download and install SIAMCAT
ENV SIAMCAT_URL "https://git.embl.de/grp-zeller/SIAMCAT/repository/development/archive.tar.gz"
RUN cd /usr/local/ \
    && wget -O siamcat.tar.gz ${SIAMCAT_URL} \
    && tar -xvzf siamcat.tar.gz \
    && mv ./SIAMCAT* ./SIAMCAT \
    && Rscript ./SIAMCAT/Rscript_flavor/00_setup.r ./siamcat.tar.gz \
    && rm siamcat.tar.gz \
    && chmod a+rx ./SIAMCAT/Rscript_flavor/*
ENV PATH="/usr/local/SIAMCAT/Rscript_flavor:${PATH}"