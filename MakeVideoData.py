import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.animation as manimation
import time
import sys
from mpl_toolkits.mplot3d import Axes3D

start_time = time.time()

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Downloads/ffmpeg-20200818-1c7e55d-win64-static/ffmpeg-20200818-1c7e55d-win64-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=30)

fig = plt.figure()

def animate(i):
    my_file = 'Adiabatic_compression_data_' + str(i) + '.txt'
    print(i)
    fig.clear()
    ax = Axes3D(fig)
    ax.set_xlim3d(-4, 4)
    ax.set_ylim3d(-4, 4)
    ax.set_zlim3d(-4, 4)
    data_x, data_y, data_z = np.genfromtxt(my_file, unpack=True)
    cont = ax.scatter(data_x, data_y, data_z, s=100, c='green')
    return cont

size_t = 2999
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('Adiabatic_compression.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
