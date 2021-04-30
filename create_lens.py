'''
Creating the facet lens mesh.

Vertices and faces are saved as a .npy file, that the built in
Python in Blender can read.

Other
-------
This was taken outside the main Blender initialization script
because the bundled in Python for Blender has no scipy-package
installed by default and installing it may be challenging.
'''

import sys
import os

import numpy as np
import scipy.spatial

import matplotlib.pyplot as plt


def create_lens(r, t, d, dxdy, close=False):
    '''
    Create a thick biconvex lens.
    
    ARGUMENTS
    r          Curvature of the both lens surfaces
    t          Thickness of the lens
    d          Lens diameter   
    
    RETURNS     vertices, faces
    '''

    hr = d/2
    
    # Hexagon side length
    ht = (2/np.sqrt(3))*hr

    x = np.linspace(-d/2, d/2, int(d/dxdy))
    y = x
    
    xv, yv = np.meshgrid(x, x)
    zv = np.nan_to_num(np.sqrt(r**2-xv**2-yv**2))
    
    vertices = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T
    
    # Circular lens capping
    vertices = vertices[vertices[:,0]**2+vertices[:,1]**2<(d/2)**2]
    
    vertices[:,2] += -r + t/2
    
    '''
    # HEXAGONAL LENS CAPPING
    inhex_vertices = []
    
    for vert in vertices:
        xi,yi,zi = vert

        # Symmetry of hexagon
        if xi < 0:
            xi = - xi
        if yi < 0:
            yi = -yi

        # Big square
        in_area1 = xi < hr and -ht/2<yi<ht/2
        # Upper triangle
        in_area2 = yi > ht/2 and yi-ht/2 < np.sqrt(ht**2-hr**2) * ((-xi/hr)+1)

        if in_area1 or in_area2:
            inhex_vertices.append([vert[0], vert[1], zi-r+t/2])

    
    inhex_vertices = np.array(inhex_vertices)
    '''
    upfaces = scipy.spatial.Delaunay(vertices[:,0:2]).simplices 
    
    faces = upfaces

    if close:
        lowfaces = scipy.spatial.Delaunay(-vertices[:,0:2]).simplices + vertices.shape[0]

        vertices = np.concatenate((vertices, -vertices), axis=0)
        faces = np.concatenate((upfaces, lowfaces), axis=0)
        
        i_corner_vertices = np.arange(0,len(vertices))[vertices[:,0] **2 + vertices[:,1] **2 > (d/2-dxdy)**2]
        
        cvertices = vertices[i_corner_vertices]
        # change to cylinder coordinates
        rho = np.sqrt(cvertices[:,0]**2+cvertices[:,1]**2)
        phi = np.arctan2(cvertices[:,1], cvertices[:,0])
        z = cvertices[:,2]
        
        scvertices = np.vstack((rho,phi,z)).T

        sidefaces = scipy.spatial.Delaunay(scvertices[:,1:]).simplices 

        global_sidefaces = []
        
        for sface in sidefaces:
            gsface = []
            for i in range(3):
                gsface.append(i_corner_vertices[sface[i]])
            global_sidefaces.append(gsface)
        
        global_sidefaces = np.array(global_sidefaces)
        
        faces = np.concatenate((faces, global_sidefaces), axis=0)

    return vertices, faces


def main():
    '''
    Create lens and save lens.npy
    '''
    
    if len(sys.argv[1:]) == 6:
        r,t,d,dxdy = [float(x) for x in sys.argv[1:5]]
        
        close = sys.argv[5]
        if close == 'True':
            close = True
        else:
            close = False
        
        savedir = sys.argv[6]
    else:
        r = 11
        t = 8
        d = 16
        dxdy = 0.3
        
        savedir = ''
        close=True
    
    os.makedirs(savedir, exist_ok=True)
    
    vertices, faces = create_lens(r, t, d, dxdy, close=close)
    np.save(os.path.join(savedir, 'lens_vertices.npy'), vertices, allow_pickle=False)
    np.save(os.path.join(savedir, 'lens_faces.npy'), faces, allow_pickle=False)


if __name__ == "__main__":
    main()
