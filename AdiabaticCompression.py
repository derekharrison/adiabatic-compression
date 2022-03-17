import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plab
import scipy.stats
from numpy import *
import matplotlib.animation as manimation
import time
import sys

#Simulation parameters
R = 0.02                           #Particle radii
L = 4.0                            #One half the length of one of the cube edges
N = 60                             #Number of particles in the system

max_vel = 9.0                      #Max value particle velocity
mass = 1.0                         #Mass of a particle
max_t = 15                         #Max simulation time
start_comp_t = 5                   #Start time compression
end_comp_t = 12                    #End time compression
v_wall = 0.2                       #Compression rate
dt = 5e-3                          #Timestep size

e  = 1.0                           #Normal restitution coefficient
mu = 0.00                          #Friction coefficient
Bo = 1.0                           #Coefficient of tangential restitution

#Initialization
err = 1e-10                        #Small number
h = L - R - err                    #Domain limit particle injection

l_wall_east = L
l_wall_west = -L
l_wall_south = -L
l_wall_north = L
l_wall_bottom = -L
l_wall_top = L

east_wall = N + 1
west_wall = N + 2
north_wall = N + 3
south_wall = N + 4
top_wall = N + 5
bottom_wall = N + 6

x = np.zeros((N))
y = np.zeros((N))
z = np.zeros((N))

I =2/5*mass*R*R
Wx = np.zeros((N))
Wy = np.zeros((N))
Wz = np.zeros((N))

avg_kin_energy_before_compression = 0
avg_kin_energy_after_compression = 0

frame_counter = 0;

#Generating initial position of first particle
x[0] = 2 * h * np.random.rand() - h
y[0] = 2 * h * np.random.rand() - h
z[0] = 2 * h * np.random.rand() - h

#Generating initial positions of remaining particles
for n in range(1,N):
    x[n] = 2 * h * np.random.rand() - h
    y[n] = 2 * h * np.random.rand() - h
    z[n] = 2 * h * np.random.rand() - h
    overlap = True;
    while overlap == True:
        overlap = False
        for i in range(0, n):
            dx = x[i] - x[n]
            dy = y[i] - y[n]
            dz = z[i] - z[n]
            if dx*dx + dy*dy + dz*dz < 1.0001*(2*R)*(2*R):
                overlap = True
        if overlap == True:
            x[n] = 2 * h * np.random.rand() - h
            y[n] = 2 * h * np.random.rand() - h
            z[n] = 2 * h * np.random.rand() - h

#Check overlap
there_is_overlap = False
for n in range(0,N):
    for i in range(n+1,N): 
        dx = x[i] - x[n]
        dy = y[i] - y[n]
        dz = z[i] - z[n]
        if dx*dx + dy*dy + dz*dz < (2*R)*(2*R):
             there_is_overlap = True
print(there_is_overlap)

#Check if in boundaries
in_boundaries = True
for n in range(0, N):
    in_boundaries = (x[n] <= l_wall_east) and (x[n] >= l_wall_west and 
                     y[n] <= l_wall_north) and (y[n] >= l_wall_south and
                     z[n] <= l_wall_top) and (z[n] >= l_wall_bottom)
    if not in_boundaries:
        in_boundaries = False

print(in_boundaries)

#Generate initial velocities
v_x = 2*max_vel*np.random.rand((N)) - np.ones((N))*max_vel
v_y = 2*max_vel*np.random.rand((N)) - np.ones((N))*max_vel
v_z = 2*max_vel*np.random.rand((N)) - np.ones((N))*max_vel

avg_kin_energy_pre = 0.0
for n in range(0,N):
    v2 = v_x[n] * v_x[n] + v_y[n] * v_y[n] + v_z[n] * v_z[n]
    avg_kin_energy_pre  = 0.5 * mass * v2 + avg_kin_energy_pre
avg_kin_energy_pre = avg_kin_energy_pre / N

#Start simulation
time = 0.0
while time <= max_t:
    print(time)
    #Calculating average kinetic energy of particle
    avg_kin_energy = 0.0
    for n in range(0,N):
        v2 = v_x[n] * v_x[n] + v_y[n] * v_y[n] + v_z[n] * v_z[n]
        avg_kin_energy  = 0.5 * mass * v2 + avg_kin_energy
    
    total_kin_energy = avg_kin_energy
    avg_kin_energy = avg_kin_energy / N

    #Save frames at regular intervals spaced at dt
    if (frame_counter == floor(time/dt)):
        my_file = 'Adiabatic_compression_data_' + str(frame_counter) + '.txt'
        file = open(my_file, "w")

        for index in range(0, size(x)):
            data_x = str(x[index])
            data_y = str(y[index])
            data_z = str(z[index])
            file.write(data_x)
            file.write(' ')
            file.write(data_y)
            file.write(' ')
            file.write(data_z)
            file.write('\n')

        file.close()
        frame_counter = frame_counter + 1
    
    if time < start_comp_t:
        avg_kin_energy_before_compression = avg_kin_energy
        kin_ratio = avg_kin_energy_pre / avg_kin_energy_before_compression
        diff = 1.0 - kin_ratio
        if abs(diff) > 1e-6:
            print('kin_ratio')
            print(kin_ratio)
            print('diff kinetic energy is too large')
            sys.exit()
    
    if time > end_comp_t:
        avg_kin_energy_after_compression = avg_kin_energy

    #Calculate minimum collision time
    collision_with_wall = False
    collision_with_particle = False
    coll_time = 1e+10

    for n in range(0,N):
        #Checking collision time between particles
        for i in range(n+1,N):
            rab = [x[n] - x[i],y[n] - y[i],z[n] - z[i]]
            vab = [v_x[n] - v_x[i],v_y[n] - v_y[i],v_z[n] - v_z[i]]
            Disc = np.dot(rab,vab)*np.dot(rab,vab)-np.dot(vab,vab)*(np.dot(rab,rab)-(2*R)*(2*R))
            coll_time_particle = dt + err

            if Disc > 0:
                coll_time_particle = ((-np.dot(rab,vab)) - sqrt(Disc)) / np.dot(vab,vab)
            else:
                coll_time_particle = 3e+8

            if coll_time_particle < coll_time and coll_time_particle >= 0:
                collision_with_particle = True
                coll_time = coll_time_particle
                coll_partner_1 = n
                coll_partner_2 = i
                
        if start_comp_t <= time and time <= end_comp_t:
            v_wall_east = -v_wall
            v_wall_west = v_wall
            v_wall_north = -v_wall
            v_wall_south = v_wall
            v_wall_top = -v_wall
            v_wall_bottom = v_wall
        else:
            v_wall_east = -0
            v_wall_west = 0
            v_wall_north = -0
            v_wall_south = 0
            v_wall_top = -0
            v_wall_bottom = 0

        coll_time_east   = (x[n] - l_wall_east + R) / (v_wall_east - v_x[n] + 1e-20)
        coll_time_west   = (x[n] - l_wall_west - R) / (v_wall_west - v_x[n] + 1e-20)
        coll_time_north  = (y[n] - l_wall_north + R) / (v_wall_north - v_y[n] + 1e-20)
        coll_time_south  = (y[n] - l_wall_south - R) / (v_wall_south - v_y[n] + 1e-20)
        coll_time_top    = (z[n] - l_wall_top + R) / (v_wall_top - v_z[n] + 1e-20)
        coll_time_bottom = (z[n] - l_wall_bottom - R) / (v_wall_bottom - v_z[n] + 1e-20)

       #Checking collision with east face        
        if coll_time_east > 0 and coll_time_east <= coll_time:
            coll_time = coll_time_east
            coll_partner_1 = n
            coll_partner_2 = east_wall
            collision_with_wall = True
            collision_with_particle = False
        
        if coll_time_west > 0 and coll_time_west <= coll_time:
            coll_time = coll_time_west
            coll_partner_1 = n
            coll_partner_2 = west_wall
            collision_with_wall = True
            collision_with_particle = False        
        
        if coll_time_north > 0 and coll_time_north <= coll_time:
            coll_time = coll_time_north
            coll_partner_1 = n
            coll_partner_2 = north_wall
            collision_with_wall = True
            collision_with_particle = False
        
        if coll_time_south > 0 and coll_time_south <= coll_time:
            coll_time = coll_time_south
            coll_partner_1 = n
            coll_partner_2 = south_wall
            collision_with_wall = True
            collision_with_particle = False        
        
        if coll_time_top > 0 and coll_time_top <= coll_time:
            coll_time = coll_time_top
            coll_partner_1 = n
            coll_partner_2 = top_wall
            collision_with_wall = True
            collision_with_particle = False         

        if coll_time_bottom > 0 and coll_time_bottom <= coll_time:
            coll_time = coll_time_bottom
            coll_partner_1 = n
            coll_partner_2 = bottom_wall
            collision_with_wall = True
            collision_with_particle = False

    if dt <= coll_time:
        x = v_x * dt * (1 - err) + x
        y = v_y * dt * (1 - err) + y
        z = v_z * dt * (1 - err) + z
        time = time + dt

        l_wall_east = v_wall_east * dt + l_wall_east
        l_wall_west = v_wall_west * dt + l_wall_west
        l_wall_north = v_wall_north * dt + l_wall_north        
        l_wall_south = v_wall_south * dt + l_wall_south
        l_wall_top = v_wall_top * dt + l_wall_top        
        l_wall_bottom = v_wall_bottom * dt + l_wall_bottom        

    elif dt > coll_time:
        x = v_x * coll_time * (1 - err) + x
        y = v_y * coll_time * (1 - err) + y
        z = v_z * coll_time * (1 - err) + z
        time = time + coll_time

        l_wall_east = v_wall_east * coll_time + l_wall_east
        l_wall_west = v_wall_west * coll_time + l_wall_west
        l_wall_north = v_wall_north * coll_time + l_wall_north        
        l_wall_south = v_wall_south * coll_time + l_wall_south
        l_wall_top = v_wall_top * coll_time + l_wall_top        
        l_wall_bottom = v_wall_bottom * coll_time + l_wall_bottom
 
        #Update velocities of colliding particles
        if collision_with_particle == True:
            ra = np.array([x[coll_partner_1],y[coll_partner_1],z[coll_partner_1]])
            rb = np.array([x[coll_partner_2],y[coll_partner_2],z[coll_partner_2]])
            va = np.array([v_x[coll_partner_1],v_y[coll_partner_1],v_z[coll_partner_1]])
            vb = np.array([v_x[coll_partner_2],v_y[coll_partner_2],v_z[coll_partner_2]])
            n_vec = (ra - rb) / sqrt(np.dot((ra-rb),(ra-rb)))
            vab=np.array([v_x[coll_partner_1]-v_x[coll_partner_2],
                          v_y[coll_partner_1]-v_y[coll_partner_2],
                          v_z[coll_partner_1]-v_z[coll_partner_2]])
            wa = np.array([Wx[coll_partner_1],
                           Wy[coll_partner_1],
                           Wz[coll_partner_1]])
            wb = np.array([Wx[coll_partner_2],
                           Wy[coll_partner_2],
                           Wz[coll_partner_2]])

            RiWi = R*wa + R*wb
            crossRiWin = np.cross(RiWi, n_vec)
            vab = vab - crossRiWin

            B1 = 7/2*(1/mass+1/mass)
            B2 = 1/mass+1/mass    

            t = (vab - n_vec*(np.dot(vab,n_vec))) / sqrt(np.dot((vab - n_vec*(np.dot(vab,n_vec))),vab - n_vec*np.dot(vab,n_vec)) + err)        
            Jn = -(1+e)*(np.dot(vab,n_vec))/B2
        
            stickyslide = (1+Bo)*np.dot(vab,t)/Jn/B1

            if mu < stickyslide: #Sliding 
                Jt = -mu*Jn              
            elif mu >= stickyslide: #Sticking
                Jt = -(1+Bo)*np.dot(vab,t)/B1      
            
            J = Jn*n_vec+Jt*t

            v_x[coll_partner_1] = J[0]/mass+v_x[coll_partner_1]
            v_y[coll_partner_1] = J[1]/mass+v_y[coll_partner_1]
            v_z[coll_partner_1] = J[2]/mass+v_z[coll_partner_1]

            v_x[coll_partner_2] = -J[0]/mass+v_x[coll_partner_2]
            v_y[coll_partner_2] = -J[1]/mass+v_y[coll_partner_2]
            v_z[coll_partner_2] = -J[2]/mass+v_z[coll_partner_2]

            crossnJ = np.cross(n_vec,-J)

            Wx[coll_partner_1] = -crossnJ[0]*R/I+Wx[coll_partner_1]
            Wy[coll_partner_1] = -crossnJ[1]*R/I+Wy[coll_partner_1]
            Wz[coll_partner_1] = -crossnJ[2]*R/I+Wz[coll_partner_1]

            Wx[coll_partner_2] = -crossnJ[0]*R/I+Wx[coll_partner_2]
            Wy[coll_partner_2] = -crossnJ[1]*R/I+Wy[coll_partner_2]      
            Wz[coll_partner_2] = -crossnJ[2]*R/I+Wz[coll_partner_2]

        elif collision_with_wall == True:
            if coll_partner_2 == east_wall:
                v_x[coll_partner_1] = 2 * v_wall_east - v_x[coll_partner_1]
            elif coll_partner_2 == west_wall:
                v_x[coll_partner_1] = 2 * v_wall_west - v_x[coll_partner_1]
            elif coll_partner_2 == north_wall:
                v_y[coll_partner_1] = 2 * v_wall_north - v_y[coll_partner_1]
            elif coll_partner_2 == south_wall:
                v_y[coll_partner_1] = 2 * v_wall_south - v_y[coll_partner_1]           
            elif coll_partner_2 == top_wall:
                v_z[coll_partner_1] = 2 * v_wall_top - v_z[coll_partner_1]
            elif coll_partner_2 == bottom_wall:
                v_z[coll_partner_1] = 2 * v_wall_bottom - v_z[coll_partner_1]


temperature_ratio_theory = (((2*L)*(2*L)*(2*L)) / ((2*l_wall_east)**3))**(2/3)
temperature_ratio_simulation = avg_kin_energy_after_compression / avg_kin_energy_before_compression

print('volume ratio V1/V2')
print(((2*L)*(2*L)*(2*L)) / ((2*l_wall_east)**3))

print('temperature_ratio_theory')
print(temperature_ratio_theory)

print('temperature_ratio_simulation')
print(temperature_ratio_simulation)
