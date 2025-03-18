import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

# Create figure with 2 3D subplots
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'scene'}]],
    subplot_titles=('Temperature Surface', 'Pressure Surface')
)

# Create a grid of x and y coordinates
x = np.linspace(-5, 5, 50)
y = np.linspace(-5, 5, 50)
X, Y = np.meshgrid(x, y)

# Create frames for animation
frames = []
num_frames = 50

for i in range(num_frames):
    # Time parameter for animation
    t = i / num_frames * 2 * np.pi
    
    # First surface: Temperature data (peaks and valleys with time evolution)
    Z1 = 3 * (1 - X/5)**2 * np.exp(-X**2/5 - (Y+1)**2/5) - \
         10 * (X/5 - X**3/5 - Y**5) * np.exp(-X**2 - Y**2) - \
         1/3 * np.exp(-(X+1)**2 - Y**2) + \
         np.sin(t) * X/2  # This part changes with time
    
    # Second surface: Pressure data (different pattern with time evolution)
    Z2 = np.sin(np.sqrt(X**2 + Y**2) + t) * np.exp(-(X**2 + Y**2)/10)
    
    # Store this frame's data
    frames.append(
        go.Frame(
            name=f'frame{i}',
            data=[
                go.Surface(
                    x=X, y=Y, z=Z1,
                    colorscale='Thermal',
                    colorbar=dict(
                        title="Temperature (°C)",
                        x=-0.07,
                        thickness=20
                    ),
                    name='temperature'
                ),
                go.Surface(
                    x=X, y=Y, z=Z2,
                    colorscale='Blues',
                    colorbar=dict(
                        title="Pressure (kPa)",
                        x=1.07,
                        thickness=20
                    ),
                    name='pressure'
                )
            ]
        )
    )

# Create initial surfaces
# First surface: Initial temperature data
Z1_init = 3 * (1 - X/5)**2 * np.exp(-X**2/5 - (Y+1)**2/5) - \
          10 * (X/5 - X**3/5 - Y**5) * np.exp(-X**2 - Y**2) - \
          1/3 * np.exp(-(X+1)**2 - Y**2)

# Second surface: Initial pressure data
Z2_init = np.sin(np.sqrt(X**2 + Y**2)) * np.exp(-(X**2 + Y**2)/10)

# Add initial traces
fig.add_trace(
    go.Surface(
        x=X, y=Y, z=Z1_init,
        colorscale='Thermal',
        colorbar=dict(
            title="Temperature (°C)",
            x=-0.07,
            thickness=20
        ),
        name='temperature'
    ),
    row=1, col=1
)

fig.add_trace(
    go.Surface(
        x=X, y=Y, z=Z2_init,
        colorscale='Blues',
        colorbar=dict(
            title="Pressure (kPa)",
            x=1.07,
            thickness=20
        ),
        name='pressure'
    ),
    row=1, col=2
)

# Update frames property
fig.frames = frames

# Configure animation settings
animation_settings = dict(
    frame=dict(duration=100, redraw=True),
    fromcurrent=True,
    transition=dict(duration=0)
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
    title='Dual Surface Plots with Synchronized Animation',
    updatemenus=updatemenus,
    sliders=sliders,
    scene=dict(
        aspectmode='cube',
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ),
    scene2=dict(
        aspectmode='cube',
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ),
    height=700,
    width=1200,
    margin=dict(l=0, r=0, b=100, t=100)
)

# Show the figure
fig.show()
