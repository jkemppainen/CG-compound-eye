'''
Algorithm(s) to create a point cloud of ommatidial locations
for the CG-compound-eye model.


'''

import os
import sys
from math import pi
import math
from math import sin, cos
import random

import numpy as np
from scipy.spatial import distance
from scipy.spatial import cKDTree as KDTree
from scipy.stats import circmean


def to_local_coordinates(x,y,z, otheta, ophi):
    '''
    Rotate so that self.cp maps to (0,1,0) in cartesian coordinates
    '''
    sin, cos = (math.sin, math.cos)
    rot = -ophi+pi/2 
    Rz = np.array([[cos(rot), -sin(rot), 0], [sin(rot), cos(rot), 0], [0,0,1]])
    rot = -pi/2+otheta
    Rx = np.array([[1,0,0], [0, cos(rot), -sin(rot)], [0, sin(rot), cos(rot)]]) 
    
    return ( Rx @ (Rz @ np.array([[x,y,z]]).T)).flatten()   


def to_global_coordinates(x,y,z, otheta, ophi):
    '''
    Rotate so that self.cp maps to (0,1,0) in cartesian coordinates
    '''
    sin, cos = (math.sin, math.cos)
    rot = ophi-pi/2 
    Rz = np.array([[cos(rot), -sin(rot), 0], [sin(rot), cos(rot), 0], [0,0,1]])
    rot = pi/2-otheta 
    Rx = np.array([[1,0,0], [0, cos(rot), -sin(rot)], [0, sin(rot), cos(rot)]])
    
    return (Rz @ Rx @ np.array([[x,y,z]]).T).flatten()


def get_xyz(R_eye, theta, phi):
    
    if callable(R_eye):
        R = R_eye(theta, phi)
    else:
        R = R_eye

    x = R * math.sin(theta) * math.cos(phi)
    y = R * math.sin(theta) * math.sin(phi)
    z = R * math.cos(theta)
    return x,y,z


def get_spherical(x,y,z):
    r = math.sqrt(x**2+y**2+z**2)
    phi = math.atan2(y, x)
    theta = math.acos(z/math.sqrt(x**2+y**2+z**2))
    return r, theta, phi


class Ommatidium:
    '''
    Hexagon having "potential" nodes.
    '''
    def __init__(self, R_eye, R_ommatidia, cp, index, parent=None):
        self.cp = np.array(cp)
        
        if callable(R_eye):
            self.R_eye = R_eye(*get_spherical(*cp)[1:])
        else:
            self.R_eye = R_eye
        
        self.R_ommatidia = R_ommatidia
        self.F = np.zeros(3, dtype=np.float)
        self.v = 0
        self.i = index
        #self.node_index

        self.vertices = 6
        #if random.random() < 0.1:
        #    self.vertices = 5

    def get_potential_nodes(self):
        nodes = []

        R, otheta, ophi = get_spherical(*self.cp)


        d_angular = 2*math.atan(self.R_ommatidia/self.R_eye)
        
        for angle in np.linspace(0, 2*pi, self.vertices+1)[:-1]:
            d_angular2 = d_angular * ( 1+(-0.5+random.random())/5 )
            
            ltheta = math.sin(angle) * d_angular2
            lphi = math.cos(angle) * d_angular2
            
            #lx,ly,lz = get_xyz(R, pi/2+ltheta, pi/2+lphi)

            ox,oy,oz = to_local_coordinates(*self.cp, otheta, ophi)
            ox += 2 * self.R_ommatidia * math.sin(angle)
            oz += 2 * self.R_ommatidia * math.cos(angle)

            nodes.append(to_global_coordinates(ox,oy,oz, otheta, ophi))

        return nodes


def initiate_star(R_eye, R_ommatidia, seed_ommatidium):
    '''
    Create starting point for create_algo algorithm.
    '''
    
    ommatidia = [seed_ommatidium]
    new_cp = (1000,1000,1000)
    
    new_ommatidia = []

    for i_direction in range(6):
        print('Creating initiation star, direction {}'.format(i_direction))
        
        new_cp = (1000,1000,1000)
        ommatidium = seed_ommatidium
        while new_cp[0] > 0:
            new_cp = ommatidium.get_potential_nodes()[i_direction]

            new_ommatidium = Ommatidium(R_eye, R_ommatidia, new_cp, len(new_ommatidia)+1)
            new_ommatidia.append(new_ommatidium)

            ommatidium = new_ommatidium

    return ommatidia+new_ommatidia


def create_algo(R_eye, R_ommatidia):
    '''
    The main eye creation algorithms
    '''
    
    initial_p = get_xyz(R_eye, pi/2, 0)
    ommatidia = [Ommatidium(R_eye(pi/2, pi/2), R_ommatidia, initial_p, 0)]
    
    ommatidia = initiate_star(R_eye, R_ommatidia, ommatidia[0])
    i_fixed = [om.i for om in ommatidia]

    gone_through = []
    go_through = None
    next_go_through = ommatidia
    
    
    while len(i_fixed) != len(ommatidia)-1:
        if go_through == next_go_through:
            break

        go_through = next_go_through
        next_go_through = []
       
        
        for ommatidium in go_through:
                        
            #if len(gone_through) > 0:
            fixed_tree = KDTree(np.array([ommatidia[ind].cp for ind in i_fixed if ind!=ommatidium.i]))

            for i_node in range(ommatidium.vertices):
                
                for i in range(10):
                    # Propose a new node point
                    new_cp = ommatidium.get_potential_nodes()[i_node]
                    #new_cp = get_xyz(R_eye, *get_spherical(*new_cp)[1:])
                    
                    distance_to_nearest = fixed_tree.query([new_cp], workers=-1)[0]

                    if distance_to_nearest < R_ommatidia:
                        break

                    if len(gone_through) == 0 or distance_to_nearest > R_ommatidia*1.25:
                        
                        new_ommatidium = Ommatidium(R_eye, R_ommatidia, new_cp, len(ommatidia))
                        ommatidia.append(new_ommatidium)
                 
                        i_fixed.append(len(ommatidia)-1)
                        next_go_through.append(new_ommatidium)
                        break

            #i_fixed.extend(neighbour_indices)
            #next_go_through.extend([ommatidia[index] for index in neighbour_indices])
            
            gone_through.append(ommatidium)
    rot = math.radians(-30)

    # Workaround, remove duplicates
    #fixed_tree = KDTree(np.array([ommatidia[ind].cp for ind in i_fixed if ind!=ommatidium.i]))
    #ommatidia = [om for om in ommatidia if fixed_tree.query([om.cp], workers=-1)[0]<R_ommatidia]

    return [(np.array([[cos(rot),0,sin(rot)],[0,1,0],[-sin(rot),0,cos(rot)]]) @ om.cp).flatten() for om in ommatidia]



def plot_cps(cps, ax=None):
    
    if ax == None:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax.clear()

    xs = [p[0] for p in cps]
    ys = [p[1] for p in cps]
    zs = [p[2] for p in cps]
    
    ax.scatter(xs, ys, zs)
    return ax



def main():

    if len(sys.argv) == 4:
        R_eye, R_ommatidia = [float(arg) for arg in sys.argv[1:-1]]
        savepath = sys.argv[3]
    elif len(sys.argv) == 5:
        R_z, R_wide, R_ommatidia = [float(arg) for arg in sys.argv[1:-1]]
        savepath = sys.argv[4]
        
        def R_eye(theta, phi):
            x,y,z = get_xyz(1, theta, phi)
            nan, theta, phi = get_spherical(x,z,y)
            return R_wide-(R_wide-R_z)*abs(math.sin(theta))
       
    else:
        #R_eye = 80s
        R_ommatidia = 8
        savepath = None
        
        R_z = 80
        R_wide = 80
        def R_eye(theta, phi):
            x,y,z = get_xyz(1, theta, phi)
            nan, theta, phi = get_spherical(x,z,y)
            return R_wide-(R_wide-R_z)*abs(math.sin(theta))

    cps = create_algo(R_eye, R_ommatidia)
    
    if not savepath:
        # Plotting
        R_max = 80
        import matplotlib.pyplot as plt
        ax = None
        ax = plot_cps(cps, ax)
        ax.set_xlim(-R_max, R_max)
        ax.set_ylim(-R_max, R_max)
        ax.set_zlim(-R_max, R_max)
        ax.set_xlabel('x') 
        ax.set_ylabel('y')        
        ax.set_zlabel('z')        
        plt.show(block=True)
        plt.pause(2)

    if savepath:
        savefn = 'curved_eye_coordinates.npy'
        np.save(os.path.join(savepath, savefn), np.array(cps), allow_pickle=False)


if __name__ == "__main__":
    main()
