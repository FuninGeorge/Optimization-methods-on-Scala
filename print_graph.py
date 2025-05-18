import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.widgets import RadioButtons
import locale

locale.setlocale(locale.LC_NUMERIC, 'de_DE.UTF-8')

def function(p):
    x, y = p
    return x ** 2 + y ** 2

def load_data(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    data = []
    for line in lines:
        parts = line.strip().split(' , ')
        x = locale.atof(parts[0])
        y = locale.atof(parts[1])
        data.append([x, y])
    return np.array(data)

data = load_data("gradient_log.csv")
x_coords, y_coords = data[:, 0], data[:, 1]
z_coords = function([x_coords, y_coords])
x = np.linspace(min(x_coords) - 0.5, max(x_coords) + 0.5, 100)
y = np.linspace(min(y_coords) - 0.5, max(y_coords) + 0.5, 100)
X, Y = np.meshgrid(x, y)
Z = function([X, Y])

fig = plt.figure(figsize=(12, 8))
ax_main = fig.add_subplot(111)

def clear_axes():
    for ax in fig.axes:
        if ax != ax_radio:
            fig.delaxes(ax)

def draw_2d():
    clear_axes()
    ax = fig.add_subplot(111)
    cs = ax.contour(X, Y, Z, levels=15, colors='black', linewidths=1)
    ax.clabel(cs, fmt="%.1f", colors='black', fontsize=10)
    ax.plot(x_coords, y_coords, '-', color='#FF4500', linewidth=2.5, alpha=0.9, label='Trajectory')
    ax.scatter(x_coords, y_coords, s=50, c='#00FFFF', edgecolor='black', linewidth=0.8, alpha=0.8, label='Temp point')
    ax.scatter(x_coords[0], y_coords[0], s=250, c='green', edgecolor='black', linewidth=2, label='Start point', zorder=3)
    ax.scatter(x_coords[-1], y_coords[-1], s=250, c='red', edgecolor='black', linewidth=2, label='End point', zorder=3)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_aspect('equal')
    ax.set_title('2D Contour', fontsize=14)
    ax.grid(alpha=0.3)
    ax.legend()
    plt.draw()

def draw_3d():
    clear_axes()
    ax3d = fig.add_subplot(111, projection='3d')
    ls = LightSource(azdeg=135, altdeg=45)
    rgb = ls.shade(Z, cmap=plt.cm.viridis, vert_exag=0.1, blend_mode='soft')
    ax3d.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=rgb, alpha=0.5)
    ax3d.plot(x_coords, y_coords, z_coords, '-', color='#FF4500', linewidth=2, label='Trajectory')
    ax3d.scatter(x_coords, y_coords, z_coords, s=30, c='#00FFFF', edgecolor='black', linewidth=1)
    ax3d.scatter(x_coords[0], y_coords[0], z_coords[0], s=200, c='green', edgecolor='black', linewidth=2)
    ax3d.scatter(x_coords[-1], y_coords[-1], z_coords[-1], s=200, c='red', edgecolor='black', linewidth=2)
    offset_z = np.min(Z) - 10
    ax3d.contour(X, Y, Z, zdir='z', offset=offset_z, levels=10, cmap='plasma', alpha=0.3)
    ax3d.set_xlim(x.min(), x.max())
    ax3d.set_ylim(y.min(), y.max())
    ax3d.set_zlim(offset_z, max(z_coords) + 5)
    ax3d.view_init(elev=35, azim=-50)
    ax3d.set_title('3D Surface', fontsize=14)
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')
    ax3d.legend()
    plt.draw()


ax_radio = plt.axes([0.02, 0.4, 0.15, 0.15])
radio = RadioButtons(ax_radio, ('2d', '3d'), active=0)

def on_select(label):
    if label == '2d':
        draw_2d()
    elif label == '3d':
        draw_3d()

radio.on_clicked(on_select)
draw_2d()
plt.show()
