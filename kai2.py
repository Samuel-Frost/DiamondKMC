import matplotlib.pyplot as plt
import numpy as np

box = [20, 20, 200]
Lx = box[0]
Ly = box[1]
time_limit = 3600 # in seconds
Nt = 10000 # number of timesteps
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

k = 2e-5 # vacancy rate constant, is the same for all of them
kNV = 2e-4 # same for kNV2

l = 2e-6 # same for V2 and V3, V4 is stable
lNV2_const = 3e-6 # NV is stable

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
print(np.sum(total))

### finite distributions
#data = []
#peepee = 10
#for n in range(Nt):
#    if n % peepee == 0:
#        data.append(V.copy())
#    V[inner] += k * laplace(V)

plot_res = Nt / 100
V_data = []
NV_data = []

for n in range(Nt):
    if n % plot_res == 0:
        V_data.append(V.copy())
        NV_data.append(V3.copy())
    dNV = kNV * (N[inner] * V[inner])
    dNV2 = kNV * (NV[inner] * V[inner])
    lNV2 = lNV2_const * (NV2[inner])

    dV2 = k * (V[inner] ** 2)
    dV3 = k * (V[inner] * V2[inner])
    lV2 = 2 * l * V2[inner] # *2 ?
    lV3 = l * V3[inner]

    V[inner] += r * laplace(V) - dV2 - dV3 - dNV - dNV2 + lV2 + lV3 
    V2[inner] += dV2 - dV3 - lV2 + lV3 
    V3[inner] += dV3 - lV3

    N[inner] -= dNV 
    NV[inner] += dNV - dNV2 + lNV2
    NV2[inner] += dNV2 - lNV2
    

total = [np.sum(V), np.sum(V2), np.sum(V3), np.sum(N), np.sum(NV), np.sum(NV2)]
#print(total)
print(np.sum(total))


### plotting
import plotly.graph_objects as go
from plotly.subplots import make_subplots

fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'scene'}]],
    subplot_titles=('NV', 'V')
)

frames = []
num_frames = len(V_data)

for i in range(len(NV_data)):
    frames.append(
        go.Frame(
            name=f'frame{i}',
            data=[
                go.Surface(
                    x=X, y=Y, z=NV_data[i],
                    colorscale='Viridis',
                    colorbar=dict(
                        title="Concentration (ppb)",
                        x=-0.07,
                        thickness=20
                    ),
                    name='temperature'
                ),
                go.Surface(
                    x=X, y=Y, z=V_data[i],
                    colorscale='Viridis',
                    colorbar=dict(
                        title="Concentration (ppb)",
                        x=1.07,
                        thickness=20
                    ),
                    name='pressure'
                )
            ]
        )
    )

fig.frames = frames


animation_settings = dict(
    frame=dict(duration=100, redraw=True),
    fromcurrent=True,
    transition=dict(duration=0)
)

fig.add_trace(
    go.Surface(
        x=X, y=Y, z=NV_data[0],
        colorscale='Viridis',
        colorbar=dict(
            title="Concentration (ppb)",
            x=-0.07,
            thickness=20
        ),
        name='temperature'
    ),
    row=1, col=1
)

fig.add_trace(
    go.Surface(
        x=X, y=Y, z=V_data[0],
        colorscale='Viridis',
        colorbar=dict(
            title="Concentration (ppb)",
            x=1.07,
            thickness=20
        ),
        name='pressure'
    ),
    row=1, col=2
)

# Configure play and pause buttons
updatemenus = [
    dict(
        type='buttons',
        showactive=False,
        buttons=[
            dict(
                label='Play',
                method='animate',
                args=[
                    None, 
                    animation_settings
                ]
            ),
            dict(
                label='Pause',
                method='animate',
                args=[
                    [None],
                    dict(
                        frame=dict(duration=0, redraw=False),
                        mode='immediate',
                        transition=dict(duration=0)
                    )
                ]
            )
        ],
        direction='left',
        pad=dict(r=10, t=10),
        x=0.1,
        y=0
    )
]

# Create and configure slider
sliders = [
    dict(
        active=0,
        yanchor='top',
        xanchor='left',
        currentvalue=dict(
            font=dict(size=16),
            prefix='Frame: ',
            visible=True,
            xanchor='right'
        ),
        transition=dict(duration=300, easing='cubic-in-out'),
        pad=dict(b=10, t=50),
        len=0.9,
        x=0.1,
        y=0,
        steps=[
            dict(
                method='animate',
                label=str(k),
                args=[
                    [f'frame{k}'],
                    dict(
                        frame=dict(duration=100, redraw=True),
                        transition=dict(duration=0),
                        mode='immediate'
                    )
                ]
            )
            for k in range(num_frames)
        ]
    )
]

# Update layout
fig.update_layout(
    title='Concentration of NV and V',
    updatemenus=updatemenus,
    sliders=sliders,
    scene=dict(
        aspectmode='cube',
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        zaxis=dict(range=[0, np.max(NV_data)])
    ),
    scene2=dict(
        aspectmode='cube',
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        zaxis=dict(range=[0, np.max(V_data)])
    ),
    height=700,
    width=1200,
    margin=dict(l=0, r=0, b=100, t=100)
)

# Show the figure
fig.show()
