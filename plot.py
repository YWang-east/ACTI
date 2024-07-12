import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = "Times New Roman"

pi = 3.14159265359

path = 'run/'

#---------------------------#
# Find the latest output    #
#---------------------------#
def latest(path):
    list_of_files = glob.glob(path + '*.csv')
    latest_file   = max(list_of_files, key=os.path.getctime)
    return latest_file

#---------------------------#
# Read case info            #
#---------------------------#
info = pd.read_csv(path + 'case_info.csv')
case = info['case'][0]
nx   = info['nx'][0]
ny   = info['ny'][0]
dim  = info['dim'][0]

#---------------------------#
# Read result data          #
#---------------------------#
data = pd.read_csv(latest(path))    # could change to read specific result file
x   = data['x'].to_numpy()
y   = data['y'].to_numpy()
rho = data['rho'].to_numpy()
u   = data['u'].to_numpy()
p   = data['p'].to_numpy()
l   = data['level'].to_numpy()

plot_quantity = rho

#---------------------------#
# Case-wise visualization   #
#---------------------------#

fig, ax = plt.subplots()

#---------------#
# mixing layer  #
#---------------#
if (case == 2):
    # For mixing layer case, extract data from diagonal (only for nx=ny)
    xd = np.zeros(nx)
    ud = np.zeros(nx)
    for i in range(nx):
        xd[i] = x[(i+1)*(ny-1)]
        ud[i] = u[(i+1)*(ny-1)] * np.sqrt(2)
    
    # analytical solution
    th = info['T_hat'][0]
    xd = 4*pi*(xd - 0.5)
    x_ex = np.linspace(-2.0*pi, 2.0*pi, 200)
    u_ex = np.exp(-th) * np.sin(x_ex)

    # plot
    ax.plot(xd, ud, 'o', markersize=3, markerfacecolor='none', c='k')
    ax.plot(x_ex, u_ex)
    ax.set_ylabel(r'$u$', fontsize=15)
    ax.grid()   
#---------------#
# 1D case       #
#---------------#
elif dim == 1:  
    ax.plot(x, plot_quantity, 'o', markersize=3, markerfacecolor='none', c='k')
    # ax.set_ylabel(r'$\rho$', fontsize=15)
#---------------#
# 2D case       #
#---------------#
elif dim == 2:
    x = np.reshape(x, (nx,ny))
    y = np.reshape(y, (nx,ny))
    r = np.reshape(plot_quantity, (nx,ny))

    plt.pcolormesh(x, y, r)
    ax.set_aspect(1)
    ax.set_xlabel(r'$x$', fontsize=15)
    ax.set_ylabel(r'$y$', fontsize=15)
    plt.colorbar()

plt.show()

