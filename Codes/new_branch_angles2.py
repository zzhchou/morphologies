#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 02:50:25 2017

@author: zzchou
"""

import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

from angle_noise_updated import angle_noise

global c_id
c_id = 1

num_morphologies = 1

inertial_weight = 3
neighbor_weight = 1
gradient_weight = 2
random_weight = 3

vertical_GF_weight = 1
radial_GF_weight = 1000
neighbor_GF_weight = 5
vertical_GF_floor = 0.9


def euclidean_distance((u1, u2, u3), (v1, v2, v3)):
    ans = math.sqrt((u1 - v1)**2 + (u2 - v2)**2 + (u3 - v3)**2)
    return ans

def cross_product((u1, u2, u3),(v1, v2, v3)):
    s1 = u2*v3 - u3*v2
    s2 = u3*v1 - u1*v3
    s3 = u1*v2 - u2*v1
    return (s1, s2, s3)


def dot_product((u1, u2, u3),(v1, v2, v3)):
    s = u1*v1 + u2*v2 + u3*v3
    return s


def magnitude_product((u1, u2, u3),(v1, v2, v3)):
    s1 = math.sqrt(u1**2 + u2**2 + u3**2)
    s2 = math.sqrt(v1**2 + v2**2 + v3**2)
    return s1*s2

def unit_vector((u1, u2, u3)):
    magnitude = euclidean_distance((u1, u2, u3), (0, 0, 0))
    
    if magnitude != 0:  
        return [u1/magnitude, u2/magnitude, u3/magnitude]
    else:
        return [0, 0, 0]

def calc_angle(vector1, vector2):
    if sum([(vector1[i] - vector2[i])**2 for i in range(0, 3)]) < 0.001:
        angle = 0
    else:
        angle = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi
    return angle

def draw_value(request):
    parameters = {

        'compartment_length': 1,
        'GF_threshold': 1,
        'random_angle': random.random()*math.pi


    }

    return parameters[request]

def plot_branches(branch_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for branch in branch_list:
        x_list.append(branch.start_coordinates[0])
        y_list.append(branch.start_coordinates[1])
        z_list.append(branch.start_coordinates[2])        
        for compartment in branch.compartment_list:
            x_list.append(compartment.end_coordinates[0])
            y_list.append(compartment.end_coordinates[1])
            z_list.append(compartment.end_coordinates[2])
        plt.plot(x_list, z_list, y_list, color='black')
        x_list = []
        y_list = []
        z_list = []

    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
#    ax.set_xlim([-100,100])
#    ax.set_ylim([-100,100])
#    ax.set_zlim([0,300])
  
    ax.set_xlim([-10,10])
    ax.set_ylim([-10,10])
    ax.set_zlim([0,20])
    
    view1 = 15
    view2 = -90
    
    ax.view_init(elev=view1, azim=view2)
    plt.show()
      
class Compartment:
    def __init__(self, start_coordinates, end_coordinates, diameter):
        
#        global c_id
        self.compartment_id = 1
        self.parent_compartment = []
#        c_id = c_id + 1
        
        self.start_coordinates = start_coordinates
        self.end_coordinates = end_coordinates
        self.diameter = diameter
        
        self.vector = [end_coordinates[i] - start_coordinates[i] for i in range(0, 3)]
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)

class Branch:
    def __init__(self, start_coordinates, start_vector, diameter):
        self.compartment_list = []
        self.start_coordinates = start_coordinates
        self.end_coordinates = start_coordinates
        self.diameter = diameter
        self.inertial_growth_vector = start_vector
        self.local_growth_vector = start_vector
        self.remote_growth_vector = start_vector
        self.x_axis = unit_vector([start_vector[1],-start_vector[0],0])
        self.canonical_end = [0,0,0]
        
        random_angle = draw_value('random_angle')
        random_vector = unit_vector([math.cos(random_angle), random.random(), math.sin(random_angle)])

        self.parent_branch = []
        self.terminal_compartment = []
        self.local_GF = 50
        
        self.grow_branch()
        
    def grow_branch(self):
        compartment_length = draw_value('compartment_length')
        
        random_angle = draw_value('random_angle')
        random_vector = unit_vector([math.cos(random_angle), random.random(), math.sin(random_angle)])        
        
        if self.terminal_compartment != []:
            self.inertial_growth_vector = self.terminal_compartment.vector
            
            
        
        self.remote_growth_vector = [self.canonical_end[i] - self.terminal_compartment.end_coordinates[i] for i in range(0, 3)]
        self.local_growth_vector = self.inertial_growth_vector
        
        remote_weight = math.sqrt(magnitude_product(self.remote_growth_vector)/magnitude_product(self.canonical_end))
        local_weight = 1 - remote_weight
        new_end_coordinates = [self.end_coordinates[i] + compartment_length*self.inertial_growth_vector[i] for i in range(0, 3)]
    
    
    
    
    
    
    
    
        new_compartment = Compartment(self.end_coordinates, new_end_coordinates, self.diameter)
        if self.compartment_list != []:
            new_compartment.parent_compartment = self.compartment_list[-1]
        elif self.parent_branch != []:
            new_compartment.parent_compartment = self.parent_branch.compartment_list[-1]
        self.compartment_list.append(new_compartment)
        self.end_coordinates = new_end_coordinates
        self.terminal_compartment = self.compartment_list[-1]
        
#        for branch2 in branch_list:
#            if branch2 != branch:
#                r = euclidean_distance(branch2.end_coordinates, branch.end_coordinates)
#                difference_vector = [branch2.terminal_compartment.vector[i] - branch.terminal_compartment.vector[i] for i in range(0, 3)]
#                
#                proj_xz = [difference_vector[0], 0, difference_vector[2]]
##                proj_xz[1] = 0
#                if proj_xz[2] >= 0:
#                    theta_list.append(math.acos(unit_vector(proj_xz)[0]))
#                else:
#                    theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))
#                    
##                theta_list.append(math.acos(difference_vector[0]/r))
#                                
#                if r != 0:
#                    total_difference_vector = [total_difference_vector[i] + difference_vector[i]/(r) for i in range(0, 3)]
#                    branch.local_GF = branch.local_GF - 300/((r**2) * (branch2.diameter**2))
#                else:
#                    branch.local_GF = 5
#                    
#                theta_list.sort()
#                                
#                biggest_angle = theta_list[0] + 2*math.pi - theta_list[-1]
#                biggest_average = (theta_list[0] + 2*math.pi + theta_list[-1])/2
#                for k in range(0, len(theta_list) - 1):
#                    current_angle = theta_list[k+1] - theta_list[k]
#                    if current_angle > biggest_angle:
#                        biggest_angle = current_angle
#                        biggest_average = (theta_list[k+1] + theta_list[k])/2
        
   
def simulate_growth(branch_list):
    for branch in branch_list:
        total_difference_vector = [0, 0, 0]
        boundary = [0, 350, 0]
        gradient_vector = unit_vector([boundary[i] - branch.end_coordinates[i] for i in range(0, 3)])
        
        
        random_angle = draw_value('random_angle')
        random_vector = unit_vector([math.cos(random_angle), random.random(), math.sin(random_angle)])
                
        r_tip = euclidean_distance([branch.end_coordinates[0], 0, branch.end_coordinates[2]], [0, 0, 0])
        
        x_max = 350
        y_max  = 375
        z_max = 200
        
        
#        branch.local_GF = 6 - 0.1*((y_boundary-branch.end_coordinates[1])/y_boundary) - 0.5*((abs(branch.end_coordinates[0]))/x_boundary) - 0.5*((abs(branch.end_coordinates[2]))/z_boundary)
        branch.local_GF = vertical_GF_weight*((y_max - branch.end_coordinates[1])/y_max - vertical_GF_floor) + radial_GF_weight*r_tip/x_max
        theta_list = []
#        print branch.end_coordinates, 1/branch.local_GF          
        
        biggest_angle = 0
        biggest_average = 0
        
        proj_xz = [branch.terminal_compartment.vector[0], 0, branch.terminal_compartment.vector[2]]
#        proj_xz[1] = 0
        if proj_xz[2] >= 0:
            theta_list.append(math.acos(unit_vector(proj_xz)[0]))
        else:
            theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))

                
        for branch2 in branch_list:
            if branch2 != branch:
                r = euclidean_distance(branch2.end_coordinates, branch.end_coordinates)
                difference_vector = [branch2.terminal_compartment.vector[i] - branch.terminal_compartment.vector[i] for i in range(0, 3)]
                
                proj_xz = [difference_vector[0], 0, difference_vector[2]]
#                proj_xz[1] = 0
                if proj_xz[2] >= 0:
                    theta_list.append(math.acos(unit_vector(proj_xz)[0]))
                else:
                    theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))
                    
#                theta_list.append(math.acos(difference_vector[0]/r))
                                
                if r != 0:
                    total_difference_vector = [total_difference_vector[i] + difference_vector[i]/(r) for i in range(0, 3)]
#                    branch.local_GF = branch.local_GF + 300/((r**2) * (branch2.diameter**2))
                    branch.local_GF = branch.local_GF + neighbor_GF_weight/(r**2)
#                    branch.local_GF = 100
                else:
                    branch.local_GF = 100
                    
                theta_list.sort()
                                
                biggest_angle = theta_list[0] + 2*math.pi - theta_list[-1]
                biggest_average = (theta_list[0] + 2*math.pi + theta_list[-1])/2
                for k in range(0, len(theta_list) - 1):
                    current_angle = theta_list[k+1] - theta_list[k]
                    if current_angle > biggest_angle:
                        biggest_angle = current_angle
                        biggest_average = (theta_list[k+1] + theta_list[k])/2
                        
                
        total_difference_vector = (math.cos(biggest_average), 0, math.sin(biggest_average))
#        total_difference_vector = unit_vector(total_difference_vector)
#        print branch.end_coordinates, 1/branch.local_GF          
        if random.random() < 1/(branch.local_GF):

#        if random.random() < float(branch.local_GF - 5)/(5*(branch.diameter**0.5)):
#            new_vector = unit_vector([inertial_weight*branch.inertial_growth_vector[i] + neighbor_weight*total_difference_vector[i] + gradient_weight*gradient_vector[i] + random_weight*random_vector[i] for i in range(0, 3)])
            new_vector = unit_vector([branch.inertial_growth_vector[i] for i in range(0,3)])
            parent_vector = unit_vector([branch.inertial_growth_vector[i] for i in range(0, 3)])
            bif_amp = 50
            
            new_vector, new_x_axis = angle_noise([0,0,0], new_vector, new_vector, branch.x_axis, 0, 0, bif_amp)
                        
            norm_vector = unit_vector(cross_product(parent_vector, [0,1,0]))
            s2 = math.sqrt(magnitude_product(new_vector, new_vector))
            s1 = s2*math.cos(bif_amp*math.pi/180)
            i_vector = [s1*parent_vector[i] for  i in range(0, 3)]
            
            p_vector = [new_vector[i] - i_vector[i] for i in range(0,3)]
            
            c_vector1 = unit_vector([p_vector[0], 0, p_vector[2]])
            c_vector2 = unit_vector([-branch.end_coordinates[0], 0, -branch.end_coordinates[2]])
            if c_vector2 == [0,0,0]:
                c_vector2 = [0,0,-1]            
            
            
            azimuth = random.gauss(30,10)
            azimuth = 0
#            c_angle = 90 - calc_angle(norm_vector, p_vector)
#            c_angle = calc_angle(c_vector1, c_vector2) + azimuth   
            c_angle = -calc_angle(norm_vector, p_vector)
            
            new_vector1, new_x_axis1 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, c_angle, 0)
                        
            new_vector2, new_x_axis2 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, -c_angle, 0)

            c_vectorA = unit_vector([new_vector1[0], 0, new_vector1[2]])
            c_vectorB = unit_vector([new_vector2[0], 0, new_vector2[2]])
            
            ang_check1 = calc_angle(c_vectorA, c_vector2)
            ang_check2 = calc_angle(c_vectorB, c_vector2)
            
            # for vertical check
            ang_check1 = new_vector2[1]
            ang_check2 = new_vector1[1]
            
#            ang_check1 = 0
#            ang_check2 = 1

#            print c_angle, c_vector2, c_vector1, c_vectorA, c_vectorB, ang_check1, ang_check2

            if ang_check1 < ang_check2:
                new_vector = new_vector1
                new_x_axis = new_x_axis1
            else:
                new_vector = new_vector2
                new_x_axis = new_x_axis2            
            
            if branch.diameter > 2:
                new_diameter = branch.diameter*0.9
                branch.diameter  = branch.diameter*0.9
            else:
                new_diameter = branch.diameter
                branch.diameter = branch.diameter
            new_branch = Branch(branch.end_coordinates, new_vector, new_diameter)
            new_branch.x_axis = unit_vector(new_x_axis)
            if branch_list != []:
                new_branch.compartment_list[-1].parent_compartment = branch.terminal_compartment
                new_branch.parent_branch = branch_list[-1]
                if new_branch.parent_branch.compartment_list != []:
                    new_branch.parent_compartment = new_branch.parent_branch.compartment_list[-1]
                else:
                    new_branch.parent_compartment = 1
            branch_list.append(new_branch)
#            print new_vector, branch.terminal_compartment.vector
            
#            if unit_vector(branch.terminal_compartment.vector)[1] > 0.95 and branch.terminal_compartment.end_coordinates[1] < 50:
#                branch.terminal_compartment.vector = [-new_vector[0], new_vector[1], -new_vector[2]]

        branch.grow_branch()
        
    return branch_list
    
# MAIN


#simulation_time = 300



origin = [0, 0, 0]
stem_diameter = 1
stem_vector = unit_vector([0, 1, 0])
branch_list = []

for j in range(0, num_morphologies):
    c_id = 1    
#    simulation_time = int(250 + random.random()*100)
    simulation_time = 20    
    
    origin = [0, 0, 0]
    stem_diameter = 1.3
    stem_vector = unit_vector([0.5, 1, 0.5])
    branch_list = []
    stem = Branch(origin, stem_vector, stem_diameter)
    branch_list.append(stem)
    
    
    for t in range(0, simulation_time):
        branch_list = simulate_growth(branch_list)
    
    
    output_title = '/home/zzchou/Documents/Data/test_output/{}{}{}'.format('morphology_output_', j+1, '.swc')
    with open(output_title, 'w') as f:
        true_origin = branch_list[0].compartment_list[0]
        line = (true_origin.compartment_id, 1, true_origin.start_coordinates[0], true_origin.start_coordinates[1], true_origin.start_coordinates[2], 4, -1)
        parent_id = 1
        strs=' '.join(str(x) for x in line)
        f.write(strs + '\n')

        for branch in branch_list:
            for compartment in branch.compartment_list:
                if compartment != []:
                    compartment.compartment_id = c_id+1
                    c_id = c_id + 1
                    if compartment.parent_compartment != []:
                        parent_id = compartment.parent_compartment.compartment_id
                    else:
                        parent_id = 1
                    line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
                    strs=' '.join(str(x) for x in line)
                    f.write(strs + '\n')
                parent_id = compartment.compartment_id
    f.closed
    
plot_branches(branch_list)




#rotation_axis = unit_vector([1,1,1])
#x_axis = unit_vector([1,-1,0])
#bif_amp = 45
#
#new_vector1, new_x_axis1 = angle_noise([0,0,0], rotation_axis, rotation_axis, x_axis, 0, 0, bif_amp)
#
#norm_vector = unit_vector(cross_product(rotation_axis, [0,1,0]))
#
#
#
#s2 = math.sqrt(magnitude_product(new_vector1, new_vector1))
#s1 = s2*math.cos(bif_amp*math.pi/180)
#i_vector = [s1*rotation_axis[i] for  i in range(0, 3)]
#
#p_vector = [new_vector1[i] - i_vector[i] for i in range(0,3)]
#
#c_angle = 90 - calc_angle(norm_vector, p_vector)
#
#print i_vector, p_vector, norm_vector, c_angle
#
#new_vector1, new_x_axis1 = angle_noise([0,0,0], new_vector1, rotation_axis, [1,0,0], 0, 0, 0)
#new_vector2, new_x_axis2 = angle_noise([0,0,0], new_vector1, rotation_axis, [1,0,0], 0, c_angle, 0)
#new_vector3, new_x_axis3 = angle_noise([0,0,0], new_vector1, rotation_axis, [1,0,0], 0, -c_angle, 0)
#
#
#
#plot1 = zip([0,0,0], new_vector1)
#plot2 = zip([0,0,0], new_vector2)
#plot3 = zip([0,0,0], new_vector3)
#plot4 = zip([0,0,0], rotation_axis)
#
#plt.figure()
#ax = plt.axes(projection='3d')
#ax.set_xlim([0,1])
#ax.set_ylim([0,1])
#ax.set_zlim([0,1])
#plt.plot(plot1[0], plot1[2], plot1[1], color='red')
#plt.plot(plot2[0], plot2[2], plot2[1], color = 'blue')
#plt.plot(plot3[0], plot3[2], plot3[1], color = 'green')
#plt.plot(plot4[0], plot4[2], plot4[1], color = 'black')