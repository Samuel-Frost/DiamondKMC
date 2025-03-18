import matplotlib.pyplot as plt
import numpy as np

box = [20, 20, 200]
Lx = box[0]
Ly = box[1]
time_limit = 3600 # in seconds
Nt = 1000 # number of timesteps
Dv = 1.5e-6 # thermal diffusivity
V0 = 1 # starting concentration in ppb 
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

k = 2e-2 # vacancy rate constant, is the same for all of them
kNV = 2e-4 # same for kNV2

l = 2e-6 # same for V2 and V3, V4 is stable
lNV2 = 3e-6 # NV is stable

### initialising distributions

V = V0 * np.exp(-(Y-10)**2/(2*1)**2) * np.exp(-X) # placeholder
# doesn't work for 0:1, only 0:2+, why?
V[:, 0:1] = 0
V_old = V.copy()
V2 = np.zeros([Nx, Ny]) # starts off as zero
V3 = np.zeros([Nx, Ny]) # starts off as zero
V4 = np.zeros([Nx, Ny]) # starts off as zero
N = 4 * np.ones([Nx, Ny])
NV = np.zeros([Nx, Ny])
NV2 = np.zeros([Nx, Ny])


### some new vectorised logic
inner = np.s_[1:-1, 1:-1] # reusable slice logic <3

def laplace(D):
    # finite difference formula vectorised
    return D[:-2, 1:-1] + D[2:, 1:-1] + D[1:-1, :-2] + D[1:-1, 2:] - 4 * D[inner]


total = [np.sum(V), np.sum(V2), np.sum(V3), np.sum(N), np.sum(NV), np.sum(NV2)]
#print(total)
#print(np.sum(total))

### finite distributions
#data = []
#peepee = 10
#for n in range(Nt):
#    if n % peepee == 0:
#        data.append(V.copy())
#    V[inner] += k * laplace(V)

plot_res = Nt / 10
V_data = []
NV_data = []

for n in range(Nt):
    if n % plot_res == 0:
        V_data.append(V.copy())
        NV_data.append(NV.copy())
    dNV = kNV * (N[inner] * V[inner])
    dNV2 = kNV * (NV[inner] * V[inner])

    dV2 = k * (V[inner] ** 2)
    dV3 = k * (V[inner] * V2[inner])

    V[inner] += k * laplace(V) - dV2 - dV3 - dNV - dNV2
    V2[inner] += dV2 - dV3
    V3[inner] += dV3

    N[inner] -= dNV 
    NV[inner] += dNV - dNV2 
    NV2[inner] += dNV2 
    

total = [np.sum(V), np.sum(V2), np.sum(V3), np.sum(N), np.sum(NV), np.sum(NV2)]
#print(total)
#print(np.sum(total))


### plotting
import plotly.graph_objects as go

fig = go.Figure(
    data=[
    go.Surface(z=NV_data[0]),
    go.Surface(z=V_data[0])
],  # Initial frame
    layout=go.Layout(
        updatemenus=[  # Add play/pause buttons
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True}],
                        "label": "Play",
                        "method": "animate",
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}],
                        "label": "Pause",
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top",
            }
        ],
        sliders=[  # Progress bar
            {
                "active": 0,
                "yanchor": "top",
                "xanchor": "left",
                "currentvalue": {"font": {"size": 20}, "prefix": "Frame: ", "visible": True, "xanchor": "right"},
                "transition": {"duration": 100, "easing": "linear"},
                "pad": {"b": 10, "t": 50},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {"args": [[str(i)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}], "label": str(i), "method": "animate"}
                    for i in range(len(V_data))
                ],
            }
        ],
    ),
    frames=[
        go.Frame(
            data=[
        go.Surface(z=NV_data[i]),
        go.Surface(z=V_data[i])
    ],
            name=str(i),
        )
        for i in range(len(V_data))
    ],
)

fig.update_layout(autosize=True)
                  #width=800, height=600)

fig.show()


quit()
#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig = plt.figure(figsize=plt.figaspect(0.5))

ax = fig.add_subplot(1, 2, 1, projection='3d')

surf1 = ax.plot_surface(X, Y, V, cmap='viridis', linewidth=0, antialiased=False)
fig.colorbar(surf1, shrink=0.5, aspect=5)

ax = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax.plot_surface(X, Y, NV, cmap='viridis', linewidth=0, antialiased=False)
surf2 = ax.plot_surface(X, Y, N, cmap='viridis', linewidth=0, antialiased=False)
fig.colorbar(surf2, shrink=0.5, aspect=5)

plt.xlabel("X")
plt.ylabel("Y")
plt.show()
