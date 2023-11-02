FROM nvidia/cuda:12.1.0-devel-ubuntu22.04

RUN apt-get update && apt-get upgrade -y

# some tools
RUN apt-get install -y git vim python3-pip

RUN apt-get install libopenmpi-dev -y 
# install python packages
RUN pip install numpy scipy cupy numba ipyparallel jupyter && \
    pip install mpi4py --force-reinstall

RUN cd /home
