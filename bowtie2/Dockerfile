FROM nathankw/centos6
#Comes with Python v2.7.10, JRE v1.8.0_91, and R v3.2.3.
#Also comes with Perl v5.10.1, which was installed as part of the "Development Tools" package.
#The directories /srv/src and /srv/software are created in the base image. 
MAINTAINER Nathaniel Watson <nathankw@stanford.edu>
RUN yum install -y bc 

#INSTALL Bowtie1.1.1
#RUN mkdir /srv/src/Bowtie1 && \
#        cd /srv/src/Bowtie1 && \
#        wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-src.zip && \
#        unzip bowtie-1.1.2-src.zip && \
#        cd bowtie-1.1.2 && \
#        make && \
#        make install
        
#INSTALL Bowtie2 2.2.8
RUN mkdir /srv/src/Bowtie2 && \
        cd /srv/src/Bowtie2 && \
        wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.8/bowtie2-2.2.8-source.zip && \
        unzip bowtie2-2.2.8-source.zip && \
        cd bowtie2-2.2.8 && \
        make && \
        make install

WORKDIR /home
ENTRYPOINT ["/usr/local/bin/bowtie2"]
LABEL version="1.0" description="Bowtie(2) docker image"

#CMD ["/bin/bash"]