# launch docker container with GPU support
docker run -it --gpus all --cap-add=SYS_ADMIN --privileged \
	-v ./:/home/ \
	-p 8888:8888 \
	hpc_python


