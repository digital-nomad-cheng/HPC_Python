# # MTF073 Computational Fluid Dynamics
# # Task 1: 2D diffusion
# # Håkan Nilsson, 2023
# # Department of Mechanics and Maritime Sciences
# # Division of Fluid Dynamics
# # Note that this is not efficient code. It is for educational purposes!

# # Clear all variables when running entire code:
# from IPython import get_ipython
# get_ipython().run_line_magic('reset', '-sf')
# # Packages needed
import numpy as np
import matplotlib.pyplot as plt

# # Close all plots when running entire code:
# plt.close('all')
# # Set default font size in plots:
plt.rcParams.update({"font.size": 12})
# import os # For saving plots

# #===================== Inputs =====================

# # Geometric and mesh inputs

L = 1.5  # Length of the domain in X direction
H = 0.5  # Length of the domain in Y direction
mI = 32 # Number of mesh points X direction.
mJ = 14 # Number of mesh points Y direction.
mesh_type = "equidistant"  # Set 'equidistant' or 'non-equidistant'
south_HN = False

# # Solver inputs

nIter = 1000  # set maximum number of iterations
resTol = 1e-3  # set convergence criteria for residuals

# #====================== Code ======================
# # Allocate arrays (np.nan used to make clear where values need to be set)
# # Note that some arrays could actually be 1D since they only have a variation
# # in one direction, but they are kept 2D so the indexing is similar for all.
nI = mI + 1  # Number of nodes in X direction, incl. boundaries
nJ = mJ + 1  # Number of nodes in Y direction, incl. boundaries
pointX = np.zeros((mI, mJ)) * np.nan  # X coords of the mesh points
pointY = np.zeros((mI, mJ)) * np.nan  # Y coords of the mesh points
nodeX = np.zeros((nI, nJ)) * np.nan  # X coords of the nodes
nodeY = np.zeros((nI, nJ)) * np.nan  # Y coords of the nodes
dx_PE = np.zeros((nI, nJ)) * np.nan  # X distance to east node
dx_WP = np.zeros((nI, nJ)) * np.nan  # X distance to west node
dy_PN = np.zeros((nI, nJ)) * np.nan  # Y distance to north node
dy_SP = np.zeros((nI, nJ)) * np.nan  # Y distance to south node
dx_we = np.zeros((nI, nJ)) * np.nan  # X size of the control volume
dy_sn = np.zeros((nI, nJ)) * np.nan  # Y size of the control volume
aE = np.zeros((nI, nJ)) * np.nan  # Array for east coefficient, in nodes
aW = np.zeros((nI, nJ)) * np.nan  # Array for wect coefficient, in nodes
aN = np.zeros((nI, nJ)) * np.nan  # Array for north coefficient, in nodes
aS = np.zeros((nI, nJ)) * np.nan  # Array for south coefficient, in nodes
aP = np.zeros((nI, nJ)) * np.nan  # Array for central coefficient, in nodes
Su = np.zeros((nI, nJ)) * np.nan  # Array for source term for temperature, in nodes
Sp = np.zeros((nI, nJ)) * np.nan  # Array for source term for temperature, in nodes
T = np.zeros((nI, nJ)) * np.nan  # Array for temperature, in nodes
k = np.zeros((nI, nJ)) * np.nan  # Array for conductivity, in nodes
k_e = np.zeros((nI, nJ)) * np.nan  # Array for conductivity at east face
k_w = np.zeros((nI, nJ)) * np.nan  # Array for conductivity at west face
k_n = np.zeros((nI, nJ)) * np.nan  # Array for conductivity at north face
k_s = np.zeros((nI, nJ)) * np.nan  # Array for conductivity at south face
res = []  # Array for appending residual each iteration
point_conv = []
F_conv = []
GE_conv = []

# # Generate mesh and compute geometric variables
if mesh_type == "equidistant":
    # # Calculate mesh point coordinates:
    for i in range(0, mI):
        for j in range(0, mJ):
            pointX[i, j] = i * L / (mI - 1)
            pointY[i, j] = j * H / (mJ - 1)

elif mesh_type == "non-equidistant":
    # This gives "symmetric log-ish around x = L/2"
    dx = np.zeros(mI)

    log_law = np.log(np.arange(int(mI/2)) + 2)
    for i in range(int(mI/2)-1):
        # If a neighbour is 15% larger we clamp it
        if log_law[i] * 1.15 < log_law[i+1]:
            log_law[i+1] = log_law[i] * 1.15
    
    dx[0:int(mI/2)] = log_law
    dx[int(mI/2)::] = np.flip(log_law)
    # We want equal spacing for the nodes in the y-direction
    dy = np.ones(mJ)
    pointX[0,:] = 0
    pointY[:,0] = 0
    for i in range(0, mI):
        for j in range(0, mJ):
            # So that the `i-1` doesn't break
            if i:
                pointX[i, j] = pointX[i-1,j] + dx[i]
            if j:
                pointY[i, j] = pointY[i,j-1] + dy[j]
    # And lastly we normalize the points
    pointX = pointX * L / pointX[-1,1]
    pointY = pointY * H / pointY[1,-1]

# # Calculate node coordinates (same for equidistant and non-equidistant):
# # Internal nodes:
for i in range(0, nI):
    for j in range(0, nJ):
        if i > 0 and i < nI - 1:
            nodeX[i, j] = 0.5 * (pointX[i, 0] + pointX[i - 1, 0])
        if j > 0 and j < nJ - 1:
            nodeY[i, j] = 0.5 * (pointY[0, j] + pointY[0, j - 1])
# # Boundary nodes:
nodeX[0, :] = 0  # Note: corner points needed for contour plot
nodeY[:, 0] = 0  # Note: corner points needed for contour plot
nodeX[-1, :] = L  # Note: corner points needed for contour plot
nodeY[:, -1] = H  # Note: corner points needed for contour plot

# # Calculate distances
# # Keep 'np.nan' where values are not needed!
for i in range(1, nI - 1):
    for j in range(1, nJ - 1):
        dx_PE[i, j] = nodeX[i + 1, j] - nodeX[i, j]
        dx_WP[i, j] = nodeX[i, j] - nodeX[i - 1, j]
        dy_PN[i, j] = nodeY[i, j + 1] - nodeY[i, j]
        dy_SP[i, j] = nodeY[i, j] - nodeY[i, j - 1]
        dx_we[i, j] = pointX[i, j] - pointX[i - 1, j]
        dy_sn[i, j] = pointY[i, j] - pointY[i, j - 1]


# # Initialize dependent variable array and Dirichlet boundary conditions
# # Note that a value is needed in all nodes for contour plot
# Our initial guess is that the temperature is 0 everywhere.
Ts, Tw = 10, 30
T = np.zeros((nI, nJ))
T[:, 0] = Ts
# If the south boundary is HN we set T to 0.
if south_HN:
    T[:, 0] = 0
T[0, :] = Tw
for j in range(0, nJ):
    T[-1, j] = 5 * (nodeY[0, j] / H - 1) + 15 * np.cos(np.pi * nodeY[0, j] / H)

# The corner values are not well-defined in the problem, they can have 2 
# different values depending on which BC you apply. Really they should be np.nan
# T[0, 0], T[-1, 0], T[0, -1], T[-1, -1] = 4 * [np.nan]

# # The following loop is for the linear solver. In the present implementation
# # it includes updates of everything that influences the linear system every
# # iteration. That may not be ideal. It may be beneficial to converge (iterate)
# # the linear solver somewhat before doing the updates. Students that are
# # interested can investigate this matter.
for iter in range(nIter):
    # # Update conductivity arrays k, k_e, k_w, k_n, k_s, according to your case:
    # # (could be moved to before iteration loop if independent of solution,
    # # but keep here if you want to easily test different cases)
    # # Keep 'np.nan' where values are not needed!
    # We need k at every cell
    kfac = 1
    for i in range(0, nI):
        for j in range(0, nJ):
            if 0.7 < nodeX[i, j] < 1.1 and 0.3 < nodeY[i, j] < 0.4:
                k[i, j] = 0.01 * kfac
            else:
                k[i,j] = 20 * kfac

    # We only need k_{s,e,n,w} at the interior nodes
    for i in range(1, nI - 1):
        for j in range(1, nJ - 1):
            fys = dy_sn[i, j] / (2 * dy_SP[i, j])
            fxe = dx_we[i, j] / (2 * dx_PE[i, j])
            fyn = dy_sn[i, j] / (2 * dy_PN[i, j])
            fxw = dx_we[i, j] / (2 * dx_WP[i, j])
            k_s[i, j] = fys * k[i, j - 1] + (1 - fys) * k[i, j]
            k_e[i, j] = fxe * k[i + 1, j] + (1 - fxe) * k[i, j]
            k_n[i, j] = fyn * k[i, j + 1] + (1 - fyn) * k[i, j]
            k_w[i, j] = fxw * k[i - 1, j] + (1 - fxw) * k[i, j]

    # # Update source term array according to your case:
    # # (could be moved to before iteration loop if independent of solution,
    # # but keep here if you want to easily test different cases)
    # # Keep 'np.nan' where values are not needed!
    # Our source is 0 everywhere, and we only need it for the interior nodes
    for i in range(1, nI - 1):
        for j in range(1, nJ - 1):
            Su[i,j] = 0

    # # Calculate coefficients:
    # # (could be moved to before iteration loop if independent of solution)
    # # Keep 'np.nan' where values are not needed!
    # # Inner node neighbour coefficients:
    # # (not caring about special treatment at boundaries):
    for i in range(1, nI - 1):
        for j in range(1, nJ - 1):
            aE[i, j] = k_e[i, j] * dy_sn[i, j] / dx_PE[i, j]
            aW[i, j] = k_w[i, j] * dy_sn[i, j] / dx_WP[i, j]
            aN[i, j] = k_n[i, j] * dx_we[i, j] / dy_PN[i, j]
            aS[i, j] = k_s[i, j] * dx_we[i, j] / dy_SP[i, j]
    b = np.zeros((nI, nJ))

    # # Modifications of aE and aW inside east and west boundaries:
    print(aW.shape)
    print(nJ-1)
    for j in range(1, nJ - 1):
        i = nI - 2  # East
        b[i,j] = aE[i,j] * T[i+1, j]
        i = 1  # West
        b[i,j] = aW[i,j] * T[i-1, j]
    # # Modifications of aN and aS inside north and south boundaries:
    # We only have modifications for the north boundary, where 
    # the difference is that we set aN to 0.
    for i in range(1, nI - 1):
        j = nJ - 2  # North
        aN[i, j] = 0
        j = 1  # South
        b[i,j] = aS[i,j] * T[i, j-1]
        # If we also wanna do homogenous neumann on the south boundary we set this to 0 as well.
        if south_HN:
            aS[i, j] = 0

    # # Inner node central coefficients:
    for i in range(1, nI - 1):
        for j in range(1, nJ - 1):
            aP[i, j] = aE[i, j] + aW[i, j] + aN[i, j] + aS[i, j]

    r = 0
    from scipy import sparse
    from scipy.sparse import linalg
    print(aP)
    print(aP[1:-1, 1:-1])
    A = sparse.diags([np.matrix.flatten(aP[1:-1, 1:-1]), -np.matrix.flatten(aN[1:-1, 1:-1])[0:-1], -np.matrix.flatten(aS[1:-1, 1:-1])[1:], -np.matrix.flatten(aE[1:-1, 1:-1]), -np.matrix.flatten(aW[1:-1, 1:-1])[nJ:]], [0, 1, -1, nJ, -nJ], format='csr')
    b = np.matrix.flatten(b[1:-1, 1:-1])
    print(A)


    # # # Solve for T using Gauss-Seidel:
    # for i in range(1, nI - 1):
    #     for j in range(1, nJ - 1):
    #         T[i, j] = (
    #             aE[i, j] * T[i + 1, j]
    #             + aW[i, j] * T[i - 1, j]
    #             + aN[i, j] * T[i, j + 1]
    #             + aS[i, j] * T[i, j - 1]
    #             + Su[i, j]
    #         ) / aP[i, j]

    
    Tsol = linalg.spsolve(A, b)
    Tsol = np.reshape(Tsol,(nI-2,nJ-2))
    T[1:-1, 1:-1] = Tsol
    break
    # # Copy T to boundaries (and corners) where homegeneous Neumann is applied:
    for i in range(1, nI - 1):
        j = nJ - 2  # North
        T[i, j + 1] = T[i, j]
        if south_HN:
            j = 0
            T[i, j] = T[i, j+1]

    # # Compute and print residuals (taking into account normalization):
    r = 0
    for i in range(1, nI - 1):
        for j in range(1, nJ - 1):
            r += abs(
                aP[i, j] * T[i, j]
                - (aE[i, j] * T[i + 1, j]
                + aW[i, j] * T[i - 1, j]
                + aN[i, j] * T[i, j + 1]
                + aS[i, j] * T[i, j - 1]
                + Su[i, j])
            )
    # Calculations for the heat flow in/out of boundaries.
    Fpos = 0
    Fneg = 0
    for j in range(1, nJ - 1):
        i = nI - 2  # East
        dT_dx = (T[i + 1, j] - T[i, j]) / dx_PE[i, j]
        qe = -k_e[i, j] * dy_sn[i, j] * dT_dx
        if qe > 0:
            Fpos += qe
        else:
            Fneg += qe
        i = 1  # West
        dT_dx = (T[i, j] - T[i - 1, j]) / dx_WP[i, j]
        qw = k_w[i, j] * dy_sn[i, j] * dT_dx
        if qw > 0:
            Fpos += qw
        else:
            Fneg += qw
    for i in range(1, nI - 1):
        # We don't calculate Fn as we know it's 0 due to the homogenous neumann condition.
        j = 1  # South
        dT_dy = (T[i, j] - T[i, j - 1]) / dy_SP[i, j]
        qs = k_s[i, j] * dx_we[i, j] * dT_dy
        if qs > 0:
            Fpos += qs
        else:
            Fneg += qs
    r /= Fpos
    
    # We only need to print every 50 iterations...
    if not iter % 50:
        print("iteration: %5d, res = %.5e" % (iter, r))
        # print(f'Fpos: {Fpos}, Fneg: {Fneg}, Fpos+Fneg: {Fpos+Fneg}, |(Fpos + Fneg) / Fp|: {abs(Fpos + Fneg) / Fpos}')
    
    # Append residual at present iteration to list of all residuals, for plotting:
    res.append(r)
    F_conv.append(Fpos)
    GE_conv.append(Fpos + Fneg)
    point_conv.append(T[-5,-5])

    # Stop iterations if converged:
    if r < resTol:
        break
print("iteration: %5d, res = %.5e" % (iter, r))
# print(f"Total heat flux over boundaries: {Fpos+Fneg:.2f} W/m, {abs(Fpos + Fneg) / Fpos * 100:.3f}% residual relative q_in")
print(f"Point convergance for x={nodeX[-5,-5]}, y={nodeY[-5,-5]}")
#================ Plotting section ================
# (only examples, more plots might be needed)

# Plot mesh
plt.figure()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Computational mesh')
plt.axis('equal')
plt.vlines(pointX[:,0],0,H,colors = 'k',linestyles = 'dashed')
plt.hlines(pointY[0,:],0,L,colors = 'k',linestyles = 'dashed')
plt.plot(nodeX, nodeY, 'ro')
plt.tight_layout()
plt.savefig('figs/mesh.pdf')

# Plot temperature contour
plt.figure()
plt.title('Temperature distribution')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.axis('equal')
tempmap=plt.contourf(nodeX.T,nodeY.T,T.T,cmap='coolwarm',levels=30)
cbar=plt.colorbar(tempmap)
cbar.set_label('Temperature [K]')
plt.tight_layout()
plt.savefig('figs/temperatureDistribution.pdf')

# Plot residual convergence
plt.figure()
plt.title('Residual convergence')
plt.xlabel('Iterations')
plt.ylabel('Residuals [-]')
resLength = np.arange(0,len(res),1)
plt.plot(resLength, res)
plt.grid()
plt.yscale('log')
plt.savefig('figs/residualConvergence.pdf')

# Plot point convergence
plt.figure()
plt.title(f'Point convergance of T in $x={nodeX[-5,-5]:.2f}, y={nodeY[-5,-5]:.2f}$')
plt.xlabel('Iterations')
plt.ylabel('T [K]')
resLength = np.arange(0,len(point_conv),1)
plt.plot(resLength, point_conv)
plt.grid()
plt.tight_layout()
plt.savefig('figs/pointConvergence.pdf')

# Plot point convergence
plt.figure()
plt.title(f'Convergence of $F$')
plt.xlabel('Iterations')
plt.ylabel('F [W/m^2]')
resLength = np.arange(0,len(F_conv),1)
plt.plot(resLength, np.abs(F_conv))
plt.grid()
plt.tight_layout()
plt.savefig('figs/FConvergence.pdf')

# Plot point convergence
plt.figure()
plt.title(f'Convergence of global energy balance')
plt.xlabel('Iterations')
plt.ylabel('$q_{sum, borders}$ [W/m^2]')
resLength = np.arange(0,len(GE_conv),1)
plt.plot(resLength, np.abs(GE_conv))
plt.grid()
plt.yscale("log")
plt.tight_layout()
plt.savefig('figs/GEConvergence.pdf')

# Plot heat flux vectors in nodes (not at boundaries)
qX = np.zeros((nI,nJ))*np.nan # Array for heat flux in x-direction, in nodes
qY = np.zeros((nI,nJ))*np.nan # Array for heat flux in y-direction, in nodes
for i in range(1, nI-1):
    for j in range(1, nJ-1):
        fys = dy_sn[i, j] / (2 * dy_SP[i, j])
        fxe = dx_we[i, j] / (2 * dx_PE[i, j])
        fyn = dy_sn[i, j] / (2 * dy_PN[i, j])
        fxw = dx_we[i, j] / (2 * dx_WP[i, j])
        T_s = fys * T[i, j - 1] + (1 - fys) * T[i, j]
        T_e = fxe * T[i + 1, j] + (1 - fxe) * T[i, j]
        T_n = fyn * T[i, j + 1] + (1 - fyn) * T[i, j]
        T_w = fxw * T[i - 1, j] + (1 - fxw) * T[i, j]
        dT_dx = (T_e - T_w) / dx_we[i, j]
        dT_dy = (T_n - T_s) / dy_sn[i, j]
        qX[i,j] = -k[i,j] * dT_dx
        qY[i,j] = -k[i, j] * dT_dy


plt.figure()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Heat flux')
plt.axis('equal')
tempmap=plt.contourf(nodeX.T,nodeY.T,T.T,cmap='coolwarm',levels=30)
plt.quiver(nodeX, nodeY, qX, qY, color="black")
cbar=plt.colorbar(tempmap)
cbar.set_label('Temperature [K]')
plt.tight_layout()
plt.savefig('figs/heatFlux.pdf')

# Plot heat flux vectors NORMAL TO WALL boundary face centers ONLY (not in corners)
# Use temperature gradient just inside domain (note difference to set heat flux)
qX = np.zeros((nI,nJ))*np.nan # Array for heat flux in x-direction, in nodes
qY = np.zeros((nI,nJ))*np.nan # Array for heat flux in y-direction, in nodes
for j in range(1,nJ-1):
    # west
    i = 0
    dT_dx = (T[i + 1, j] - T[i, j]) / dx_WP[i + 1,j]
    qX[i,j] = -k[i,j] * dT_dx
    qY[i,j] = 0
    # east
    i = nI-1
    dT_dx = (T[i, j] - T[i - 1, j]) / dx_PE[i - 1,j]
    qX[i,j] = -k[i,j] * dT_dx
    qY[i,j] = 0
for i in range(1,nI-1):
    # south
    j = 0
    dT_dy = (T[i, j + 1] - T[i, j]) / dy_SP[i, j + 1]
    qX[i,j] = 0
    qY[i,j] = -k[i,j] * dT_dy
    # north
    j = nJ-1
    dT_dy = (T[i, j] - T[i, j - 1]) / dy_SP[i, j - 1]
    qX[i,j] = 0
    qY[i,j] = -k[i,j] * dT_dy
plt.figure()
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Wall-normal heat flux \n (from internal temperature gradient)')
plt.axis('equal')
tempmap=plt.contourf(nodeX.T,nodeY.T,T.T,cmap='coolwarm',levels=30)
cbar=plt.colorbar(tempmap)
cbar.set_label('Temperature [K]')
plt.quiver(nodeX, nodeY, qX, qY, color="black")
plt.xlim(-0.5*L, 3/2*L)
plt.ylim(-0.5*H, 3/2*H)
plt.tight_layout()
plt.savefig('figs/wallHeatFlux.pdf')
