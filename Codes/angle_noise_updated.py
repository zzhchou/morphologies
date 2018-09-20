from __future__ import division

import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def euclidean_distance((u1, u2, u3), (v1, v2, v3)):
    ans = math.sqrt((u1 - v1)**2 + (u2 - v2)**2 + (u3 - v3)**2)
    return ans

def cross_product((u1, u2, u3),(v1, v2, v3)):
    s1 = u2*v3 - u3*v2
    s2 = u3*v1 - u1*v3
    s3 = u1*v2 - u2*v1
    return (s1, s2, s3)

def angle_noise(origin, end_point, y_axis, x_axis, x_rot, y_rot, z_rot):
    
    origin_array = np.asarray(origin)
    end_array = np.asarray(end_point)
    x_axis_array = np.asarray(x_axis)
    
    vector_real = end_array - origin_array
    vector = np.asarray(y_axis)
    mag = np.sqrt(vector.dot(vector))
    
    vector_unit = vector/mag
    vector_list = np.ndarray.tolist(vector_unit)
    
    ux = vector_list[0]
    uy = vector_list[1]
    uz = vector_list[2]
    
    perp_vector_list = x_axis
    
    vx = perp_vector_list[0]
    vy = perp_vector_list[1]
    vz = perp_vector_list[2]
    
    perp_vector_list2 = cross_product(vector_list, perp_vector_list)
    
    wx = perp_vector_list2[0]
    wy = perp_vector_list2[1]
    wz = perp_vector_list2[2]


    rot_x = x_rot*math.pi/180
    rot_y = y_rot*math.pi/180
    rot_z = z_rot*math.pi/180


#    print rot_x, rot_y, rot_z

#    s = math.sin(rot_x)
#    c = math.cos(rot_x)
#    omc = 1-c

    R_x = np.matrix(( (math.cos(rot_x) + vx*vx*(1-math.cos(rot_x)), vx*vy*(1-math.cos(rot_x))-vz*math.sin(rot_x), vx*vz*(1-math.cos(rot_x))+vy*math.sin(rot_x) ),
                      (vy*vx*(1-math.cos(rot_x))+vz*math.sin(rot_x), math.cos(rot_x)+vy*vy*(1-math.cos(rot_x)), vy*vz*(1-math.cos(rot_x))-vx*math.sin(rot_x)),
                      (vz*vx*(1-math.cos(rot_x))-vy*math.sin(rot_x), vz*vy*(1-math.cos(rot_x))+vx*math.sin(rot_x), math.cos(rot_x)+vz*vz*(1-math.cos(rot_x))) ))
    
    
    R_y = np.matrix( ((math.cos(rot_y) + ux*ux*(1-math.cos(rot_y)), ux*uy*(1-math.cos(rot_y))-uz*math.sin(rot_y), ux*uz*(1-math.cos(rot_y))+uy*math.sin(rot_y)),
                      (uy*ux*(1-math.cos(rot_y))+uz*math.sin(rot_y), math.cos(rot_y)+uy*uy*(1-math.cos(rot_y)), uy*uz*(1-math.cos(rot_y))-ux*math.sin(rot_y)),
                      (uz*ux*(1-math.cos(rot_y))-uy*math.sin(rot_y), uz*uy*(1-math.cos(rot_y))+ux*math.sin(rot_y), math.cos(rot_y)+uz*uz*(1-math.cos(rot_y)))) )
   
    
    R_z = np.matrix( ((math.cos(rot_z) + wx*wx*(1-math.cos(rot_z)), wx*wy*(1-math.cos(rot_z))-wz*math.sin(rot_z), wx*wz*(1-math.cos(rot_z))+wy*math.sin(rot_z)),
                      (wy*wx*(1-math.cos(rot_z))+wz*math.sin(rot_z), math.cos(rot_z)+wy*wy*(1-math.cos(rot_z)), wy*wz*(1-math.cos(rot_z))-wx*math.sin(rot_z)),
                      (wz*wx*(1-math.cos(rot_z))-wy*math.sin(rot_z), wz*wy*(1-math.cos(rot_z))+wx*math.sin(rot_z), math.cos(rot_z)+wz*wz*(1-math.cos(rot_z)))) )

#    print R_x
#    print R_y
#    print R_z

    R = R_x * R_y * R_z
    
    v_array = np.asarray(vector_real)

    new_vector = np.ndarray.tolist(np.dot(v_array, R))[0]
    new_x_axis = np.ndarray.tolist(np.dot(x_axis_array, R))[0]
    
    new_end = [origin[i] + new_vector[i] for i in range(0, 3)]
    return new_end, new_x_axis
    
# MAIN
    
#origin = [0,0,0]
#vector = [0.5, 1, 0.5]
#
#theta = 30
#phi = 20
#new_vector = angle_noise(origin, vector, theta, phi)
#
#vector_coordinates = zip(origin, vector)
#new_vector_coordinates = zip(origin, new_vector)
#
#x1 = vector_coordinates[0]
#y1 = vector_coordinates[1]
#z1 = vector_coordinates[2]
#
#x2 = new_vector_coordinates[0]
#y2 = new_vector_coordinates[1]
#z2 = new_vector_coordinates[2]
#
#
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#
#plt.plot(x1, z1, y1)
#plt.plot(x2, z2, y2)
#ax.set_xlim([-1,1])
#ax.set_ylim([-1,1])
#ax.set_zlim([0,2])
#plt.show()