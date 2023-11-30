import logging

import numpy as np

import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)

def generate_input_data(grid_size: int):
    logger.info("generate input data with grid size: {}".format(grid_size))
    ni=grid_size
    nj=grid_size
    yfac=1.05 # stretching
    ymax=2
    xmax=3
    viscos=1/1000
    dy=0.1
    yc=np.zeros(nj+1)
    yc[0]=0.
    for j in range(1,int(nj/2)+1):
        yc[j]=yc[j-1]+dy
        dy=yfac*dy


    ymax_scale=yc[int(nj/2)]

    # cell faces
    for j in range(1,int(nj/2)+1):
       yc[j]=yc[j]/ymax_scale
       yc[nj-j+1]=ymax-yc[j-1]

    yc[int(nj/2)]=1


    logger.info('y+ {}'.format(0.5*yc[1]/viscos))
    # make it 2D
    y2d=np.repeat(yc[None,:], repeats=ni+1, axis=0)

    y2d=np.append(y2d,nj)
    np.savetxt('y2d.dat', y2d)

    # x grid
    xc=yc
    # make it 2D
    x2d=np.repeat(xc[:,None], repeats=nj+1, axis=1)
    x2d_org=x2d
    x2d=np.append(x2d,ni)
    np.savetxt('x2d.dat', x2d)

    # check it
    datay= np.loadtxt("y2d.dat")
    y=datay[0:-1]
    nj=int(datay[-1])

    y2=np.zeros((ni+1,nj+1))
    y2=np.reshape(y,(ni+1,nj+1))

    datax= np.loadtxt("x2d.dat")
    x=datax[0:-1]
    ni=int(datax[-1])

    x2=np.zeros((ni+1,nj+1))
    x2=np.reshape(x,(ni+1,nj+1))

