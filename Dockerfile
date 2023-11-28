FROM nvidia/cuda:12.1.0-devel-ubuntu22.04

RUN apt-get update && apt-get upgrade -y

# some tools
RUN apt-get install -y git vim python3-pip

RUN apt-get install libopenmpi-dev -y 
# install python packages
RUN pip install numpy scipy cupy numba ipyparallel jupyter matplotlib pandas && \
    pip install mpi4py --force-reinstall
RUN pip install https://github.com/mlc-ai/package/releases/download/v0.9.dev0/mlc_ai_nightly_cu121-0.12.dev1880-cp310-cp310-manylinux_2_28_x86_64.whl
RUN cd /home
