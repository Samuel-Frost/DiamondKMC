import matplotlib.pyplot as plt
import numpy as np

box = [20, 20, 200]
Lx = box[0]
Ly = box[1]
time_limit = 3600 # in seconds
Nt = 500 # number of timesteps
Dv = 1.5e-6 # thermal diffusivity
V0 = 100 # starting concentration in ppb 
alpha = 1000 # raman scattering coeff
wZ = 0.404272342 # in micrometres

# same density throughout
density = 5
Nx = box[0] * density 
Ny = box[1] * density
Nz = box[2] * density
x = np.linspace(0, box[0], Nx)
y = np.linspace(0, box[1], Ny)
z = np.linspace(0, box[2], Nz)

dx = Lx / (Nx) 
dy = Ly / (Ny) 
dt = time_limit / Nt

r = Dv * dt / dy**2
if r > 0.5:
    raise Exception("Solution is unstable") 

[X, Y, Z] = np.meshgrid(x, y, z) # 3D array
[X, Y] = np.meshgrid(x, y) # 3D array

k = 2e-5 # vacancy rate constant
kNV = 2e-4

### initialising distributions

V = V0 * np.exp(-(Y-10)**2/(2*2)**2) * np.exp(-X) # placeholder
# doesn't work for 0:1, only 0:2+, why?
V[:, 0:2] = 0
V_old = V.copy()
V2 = np.zeros([Nx, Ny]) # starts off as zero
N = 4 * np.ones([Nx, Ny])
NV = np.zeros([Nx, Ny])

print(np.sum(V), np.sum(NV), np.sum(N))
### finite distributions
for n in range(Nt):
    V_new = V
    NV_new = NV
    N_new = N
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            # is k here but maybe should be r?
            V_new[i, j] = V[i, j] + k * (V[i+1, j] + V[i, j+1] + V[i-1, j] + V[i, j-1] - 4 * V[i, j]) - kNV*(N[i, j]*V[i, j])
            NV_new[i, j] = NV[i, j] + kNV*(N[i, j] * V[i, j])
            N_new[i, j] = N[i, j] - kNV*(N[i, j] * V[i, j])
    V = V_new 
    NV = NV_new
    N = N_new



print(np.sum(V), np.sum(NV), np.sum(N))
### plotting

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig = plt.figure(figsize=plt.figaspect(0.5))

ax = fig.add_subplot(1, 2, 1, projection='3d')

surf1 = ax.plot_surface(X, Y, V, cmap='viridis', linewidth=0, antialiased=True)
fig.colorbar(surf1, shrink=0.5, aspect=5)

ax = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax.plot_surface(X, Y, NV, cmap='viridis', linewidth=0, antialiased=True)
surf2 = ax.plot_surface(X, Y, N, cmap='viridis', linewidth=0, antialiased=True)
fig.colorbar(surf2, shrink=0.5, aspect=5)

plt.xlabel("X")
plt.ylabel("Y")
plt.show()
