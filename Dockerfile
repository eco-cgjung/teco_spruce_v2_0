FROM centos:7
MAINTAINER Mark Stacy <markstacy@ou.edu>

RUN mkdir /source
COPY . /source
WORKDIR /source
RUN yum install -y gcc gcc-gfortran
RUN gfortran -o /source/teco_spruce_v2_0 /source/TECO_SPRUCE_2_0.f90
RUN chmod +x /source/teco_spruce_v2_0
ENTRYPOINT ["./teco_spruce_v2_0"]
CMD ["./input/SPRUCE_pars.txt", "./input/SPRUCE_forcing.txt", "./input/SPRUCE_obs.txt", "/source/output/", "0"]
