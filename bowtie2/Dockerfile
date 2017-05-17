############################################################
# Dockerfile to build Bowtie 2 container images
# Based on Ubuntu
############################################################
# Set the base image to Ubuntu
FROM ubuntu:16.04
# File Author / Maintainer
MAINTAINER Derrick DeConti <deconti@jimmy.harvard.edu>
# Update the repository sources list
RUN apt-get update
# Install compiler and perl stuff
RUN apt-get install --yes build-essential gcc-multilib apt-utils zlib1g-dev git
# Download and build 
RUN mkdir /srv/src/Bowtie2 && \
        cd /srv/src/Bowtie2 && \
        wget https://github.com/BenLangmead/bowtie2/archive/v2.3.2.zip && \
        unzip bowtie2-2.3.2.zip && \
        cd bowtie2-2.3.2 && \
        make && \
        make install
ENTRYPOINT []
LABEL version="2.3.2" description="Bowtie2"
CMD ["/bin/bash"]