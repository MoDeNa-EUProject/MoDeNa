FROM ubuntu:16.04

WORKDIR /app
RUN apt-get update && apt-get install -y automake libltdl-dev libltdl7 mongodb  python-pip python-scipy python-rpy2 python-blessings r-base r-base-dev build-essential swig cmake

COPY . .
RUN pip install -r requirements.txt
RUN cd src && \
  cmake -DCMAKE_INSTALL_PREFIX:PATH=/app/MoDeNa . \
  && make \
  && make install
  
