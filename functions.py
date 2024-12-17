import numpy as np
from scipy import interpolate

def getnearpos(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def parametrization_axis_ratio(diameters):

    d_stop1, d_stop2 = 1.35, 4.4  # mm

    ax1 = np.array(
        [1 / (0.9939 + 0.00736 * D - 0.018485 * D ** 2 + 0.001456 * D ** 3) for D in diameters])  # keenan et. al 2001
    ax2 = np.array([1 / (1.0048 + 0.0057 * D - 2.628 * D ** 2 + 3.682 * D ** 3 - 1.677 * D ** 4) for D in
                            diameters / 10])  # Bringi and Chandrasekar 2001

    x = np.array((4, 5, 6, 7, 8))
    y = np.array((0.778, 0.708, 0.642, 0.581, 0.521))
    f = interpolate.interp1d(x, 1/y, fill_value="extrapolate") # (x,y)
    ax3 = f(diameters) # Beard and Chuang 1987

    axis_ratio = np.concatenate((ax1[:getnearpos(diameters,d_stop1)],
                         ax2[getnearpos(diameters,d_stop1):getnearpos(diameters,d_stop2)],
                         ax3[getnearpos(diameters,d_stop2):]))

    return axis_ratio


def velocity_diameters_rel(d_m):  # d in meters
    global c
    c = 1.2 * 1e8

    Vcloud = c * np.power(d_m / 2, 2)
    Vdrizzle = 8333 * d_m / 2 - 0.0833
    Vrain = 9.65 - 10.3 * np.exp(-.6e3 * d_m)

    def heavirange(x, x1, x2):
        return np.multiply(np.heaviside(x - x1, 0.5), (1 - np.heaviside(x - x2, 0.5)))

    V1 = np.multiply(Vcloud, heavirange(d_m, 0, 2 * 57.33 * 1e-6))

    V2 = np.multiply(Vdrizzle, heavirange(d_m, 2 * 57.33 * 1e-6, 2 * 430 * 1e-6))

    V3 = np.multiply(Vrain, heavirange(d_m, 2 * 430 * 1e-6, 9))

    return V1 + V2 + V3

import os
import matplotlib.pyplot as plt


def fast_plotting(x, y, color, xlabel="X-axis", ylabel="Y-axis", title="Plot Title", save_dir="./", filename="plot.png", ylimits=None):
    """
    Creates a quick plot with the given x and y values and saves it to the specified directory.

    Parameters:
    x (list or array-like): Data for the x-axis.
    y (list or array-like): Data for the y-axis.
    xlabel (str): Label for the x-axis. Default is "X-axis".
    ylabel (str): Label for the y-axis. Default is "Y-axis".
    title (str): Title of the plot. Default is "Plot Title".
    save_dir (str): Directory where the plot will be saved. Default is the current directory.
    filename (str): Name of the file to save the plot. Default is "plot.png".

    Returns:
    None
    """
    ticksize = 14
    labelsize = 18


    plt.figure()
    plt.plot(x, y, color=color)  # Customize the plot style as needed
    plt.xlabel(xlabel, size= labelsize)
    plt.ylabel(ylabel, size= labelsize)
    plt.title(title, size= labelsize-2)
    plt.grid(ls='--', lw=0.5)

    plt.xlim(0,10)
    if ylimits:
        plt.ylim(ylimits)

    plt.xticks(size=ticksize)
    plt.yticks(size=ticksize)

    save_path = os.path.join(save_dir, filename)
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight', dpi=200)
    plt.show()
    plt.close()

    print(f"Plot saved at: {save_path}")
