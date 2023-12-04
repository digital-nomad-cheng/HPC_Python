import numpy as np
import sys


args = sys.argv
if len(args) > 1:
    scalefactor = int(args[1])
else:
    scalefactor = 1

ni = 200 * scalefactor
nj = 100 * scalefactor


xfac = 1.005
xfac = 1.01
viscos = 3.57e-5
dy = 7.83208e-004
yc = np.zeros(nj + 1)
yc[0] = 0.0
yfac = 1.1
ustar = 1 / 25
for j in range(1, nj + 1):
    yc[j] = yc[j - 1] + dy
    yplus = yc[j] * ustar / viscos
    if dy < 0.05:
        dy = yfac * dy
    # print("j=%d, y=%7.2E, yplus=%7.2E, dy=%7.2E" % (j, yc[j], yplus, dy))

ymax_scale = yc[nj]

# u_inf=1 => ustar=1.25
ustar = 1 / 25

# make it 2D
y2d = np.repeat(yc[None, :], repeats=ni + 1, axis=0)

y2d = np.append(y2d, nj)
np.savetxt("y2d.dat", y2d)

#
xc = np.zeros(ni + 1)
dx = 0.03
for i in range(1, ni + 1):
    xc[i] = xc[i - 1] + dx
    if dx < 0.5:
        dx = dx * xfac
    # print("i=%d, x=%7.2E, dx=%7.2E" % (i, xc[i], dx))


# make it 2D
x2d = np.repeat(xc[:, None], repeats=nj + 1, axis=1)
x2d_org = x2d
x2d = np.append(x2d, ni)
np.savetxt("x2d.dat", x2d)


# check it
datay = np.loadtxt("y2d.dat")
y = datay[0:-1]
nj = int(datay[-1])

y2 = np.zeros((ni + 1, nj + 1))
y2 = np.reshape(y, (ni + 1, nj + 1))

datax = np.loadtxt("x2d.dat")
x = datax[0:-1]
ni = int(datax[-1])

x2 = np.zeros((ni + 1, nj + 1))
x2 = np.reshape(x, (ni + 1, nj + 1))
