''' Realistic compound eye CG model with optics and moving rhabdomeres.


This script opened in Blender's Text Editor can be ran by pressing the Run Script button,
leading to a creation of the CG eye model.

The initialization part is slow (see TODO below).


Dependencies
-------------
    - Blender >= 2.8 (tested with 2.8.1)
    - Luxcorerender >= 2.3 (tested with 2.3)
    - Working Python installation with Scipy


Parameters
----------
    The CG model has easily modifiable parameters to give different outcomes
    
    IMPORTANT PARAMETERS
    --------------------
        EXTERNAL_PYTHON : string
            In path name of the external python exetable, for example 
            "python3.7" or "/usr/bin/python" etc.
                Default : "python"
        EYE_DESIGN : string
            Name one of the EYE_PARAMETER classes, 'drosophila_eye' or 'generic_eye'
            available atm.
                Default : "drosophila_eye"
      

    EYE PARAMETERS
    --------------
        For details, see generic_eye -class in the code below.
        
        See also the main function arguments for extra options.


    PERFORMANCE TUNING PARAMETERS (lines 50-75)
    -------------------------------------------
        LENS_DXDY : float
            Amount of detail on the lens surfaces. Higher values lead
            to less details and lower to more.
                Default: 0.3
        ONLY_RHABDOMERES : bool
            If True, skip generating the lenses and the screening pigments
                Default: False
        SCALE_FACTOR : float
            Model parameters are specified in real life micrometers, but by default
            generater models in Blender meter units. For the whole eye, this would
            lead to very large eye by default, and rendering very large things
            can lead to artifacts.
            Smaller values make the generated structures smaller.
                Default : 0.001 (for single ommatidium shots try 1.)
        RHABDOMERE_VERTICES : int
            Number of corners for the rhabdomere cylinders. Higher values
            make rhabdomeres cross-section more circular.
                Default : 16
        LOWLEVEL_OPS : bool
            (experimental) If True, avoid bpy.ops and link as late as possible,
            shortening the model init time from tens of minutes -> under a minute.
            Animations do not work yet with this and rhabdomeres not properly rotated.
    

TODO
-----
    - Create materials automatically beforehand
    - Make the model initialiaztion faster (avoiding bpy.ops in init)
    - Error handling with the external python calls or remove external
        and install scipy in Blender?
    - Convert this script into a Blender Plugin?


License
--------
This script is licenced under the GPL-3.0 License.

For the copyrigh of the generated output (the 3D model), please refer to
https://www.gnu.org/licenses/gpl-faq.en.html#GPLOutput
'''

import bpy
import os
import json
import subprocess
from math import radians, pi, sin, cos, atan, atan2, sqrt, asin, acos, tan
import copy
import random
import numpy as np
import mathutils
import time



# PERFORMANCE TUNING PARAMETERS
# -----------------------------
LENS_DXDY = 0.3
ONLY_RHABDOMERES = False        # If True, omit screening pigments and the lenses
SCALE_FACTOR = 0.001            # Scaling factor from parameter micrometers to Blender's units
RHABDOMERE_VERTICES = 16        # Amount of circularity in rhabdomere cylinders
LOWLEVEL_OPS = True             # Avoid bpy.ops and link as late as possible

# OTHER
# --------
EXTERNAL_PYTHON="python"        # In path name of the external python exetable
EYE_DESIGN = 'drosophila_eye'   # 'drosophila_eye' or 'generic_eye' available atm.


# EYE PARAMETERS (using classes for namespacing)
# ----------------------------------------------
class generic_eye:
    '''
    Parameters for a generic eye (faintly resembling the eye of Drosophila, less ommatidia)
    for generic simulations / tests.
    
    The following list contains complete description of parameters avaiable (and having an
    effect) by default.
    
    ------------------------------------------------------------------------------------------
    PARAMETER               DESCRIPTION                                         UNITS (if any)
    ------------------------------------------------------------------------------------------
    
    OBJECT_PLANE_TYPE       "rhabdomeres"n
    
    # (rhabdomeres)
    RHABDOMERE_LOCATIONS    A list of (x,y,z)-coordinates                       um
    RHABDOMERE_DIAMETERS    Rhabdomere diameters,                               um
    RHABDOMERE_ROTATION     Rotation along z-axis                               radians
    RHABDOMERE_DEPTH        Lenght of the rhabdomeres in z                      um
    
    # (screening pigments)
    HEXAGON_VERTICES
    HEXAGON_DIMENSIONS
    HEXAGON_ROTATION
    HEXAGON_INOUT_RATIO
    
    # (lens parameters)
    LENS_D                  Lens diameter                                       um
    LENS_R                  Curvature of the lens surfaces (both same)          um
    LENS_T                  Lens thickness                                      um
    IOR_OUTER
    IOR_INNER
    LENS_TO_TIPS_DISTANCE   Distance from the lens CP to the rhabdomere tips    microns
    
    # (compound eye)
    CURVED_EYE_METHOD       "icosphere" or "external"
    EYE_RADIUS              The main radius of the eye                          um
    EYE_HORIZONTAL_R        Eye radius in (x,y) plane                           um
    OMMATIDIA_LIMIT         Stop creating ommatidia after this (failsafe)
    EYE_LOWER_ANGLE         Maximum OA angle in coronal plane from top          radians
    
    # (two eyes / multieye)
    EYE_HEAD_CP_DISTANCE    Distance between the L and R eye inner corners      um
    EYE_EXTEND_INVARDS      Amount of eyes extending inwards or outwards        
    
    ------------------------------------------------------------------------------------------
    ------------------------------------------------------------------------------------------
    '''
    OBJECT_PLANE_TYPE = "rhabdomeres"

    RHABDOMERE_LOCATIONS = [(-1.6881, 1.0273, 0.0000), (-1.8046, -0.9934, 0.0000),
        (-1.7111, -2.9717, 0.0000), (-0.0025, -1.9261, 0.0000),
        (1.6690, -0.9493, 0.0000), (1.6567, 0.9762, 0.0000),
        (0.0045, -0.0113, 0.0000)]
    RHABDOMERE_DIAMETERS = [1.8627,1.8627,1.8627,1.8627,1.8627,1.8627, 1.5743]
    RHABDOMERE_ROTATION = radians(-45)
    RHABDOMERE_DEPTH = 3
    
    HEXAGON_VERTICES = 6
    HEXAGON_DIMENSIONS = (15.99, 16.477, 1)
    HEXAGON_ROTATION = radians(0)
    HEXAGON_INOUT_RATIO = 0.95 

    LENS_D = 16
    LENS_R = 11
    LENS_T = 8
    IOR_OUTER = 1.45
    IOR_INNER = 1.45/1.34
    LENS_TO_TIPS_DISTANCE = 27
        
    # Parameters of the eye
    CURVED_EYE_METHOD = "icosphere" 
    EYE_RADIUS =  LENS_D/2*10
    EYE_HORIZONTAL_R = None
    EYE_SCALING = [1,1,1]
    OMMATIDIA_LIMIT=100

    EYE_LOWER_ANGLE = None
    MIRROR_HEMISPHERES = False
    EYE_HEAD_CP_DISTANCE = 0



class drosophila_eye(generic_eye):
    '''
    Additional overriding parameters to generate Drosophila's eye more realisticly.
    
    See generic_eye class for documentation and base values.
    '''
    
    # Rhabdomere diameters from their areas, digitized from 
    # Juusola et al, Microsaccadic sampling (2017), Appendix 5 fig 1
    RHABDOMERE_DIAMETERS = [2*sqrt(A/pi) for A in [2.6797385620915035, 2.0980392156862746,
        2.2287581699346406, 1.9411764705882355, 2.1764705882352944, 2.65359477124183]]
    
    # Based on ranceschini Pupil and Pseudopupil in the Compound Eye of Drosophila (1972),
    # scale R7/8 down from R1-6 average, relative
    RHABDOMERE_DIAMETERS.append(np.mean(RHABDOMERE_DIAMETERS)*(1.5743/1.8627))
    
    RHABDOMERE_ROTATION = radians(0)
    HEXAGON_DIMENSIONS = [16, 16, 1]
    LENS_TO_TIPS_DISTANCE = 21
    LENS_D = 16 
    
    EYE_RADIUS =  1.1*1000*(400*(0.8/985))/2
    EYE_HORIZONTAL_R = 1.1*1000*(470*(0.8/985))/2
    
    OMMATIDIA_LIMIT=800
    CURVED_EYE_METHOD = "external" 
    EYE_LOWER_ANGLE = radians(120/2)
    MIRROR_HEMISPHERES = True
    EYE_HEAD_CP_DISTANCE = 1000*(476*(0.8/985))/2
    EYE_EXTEND_INVARDS = 0




def rescale1(factor):
    '''
    Scale all the params by the factor. 
    '''
    global LENS_DXDY
    LENS_DXDY *= factor
    
    params.RHABDOMERE_LOCATIONS = np.array(params.RHABDOMERE_LOCATIONS) * factor
    params.RHABDOMERE_DIAMETERS = np.array(params.RHABDOMERE_DIAMETERS) * factor
    params.HEXAGON_DIMENSIONS = np.array(params.HEXAGON_DIMENSIONS ) * factor
    params.HEXAGON_DIMENSIONS[2] = 1

    params.LENS_D *= factor
    params.LENS_R *= factor
    params.LENS_T *= factor
    params.LENS_TO_TIPS_DISTANCE *= factor
    params.EYE_RADIUS *= factor
    params.EYE_HEAD_CP_DISTANCE *= factor
    params.RHABDOMERE_DEPTH *= factor



def set_material(object, material_name, append=False):
    '''
    Setting material to a Blender object.
    '''
    
    material = bpy.data.materials.get(material_name)
    if material == None:
        return None

    if append:
        object.data.materials.append(material)
        return None

    # Set material to the object
    if object.data.materials:
        object.data.materials[0] = material
    else:
        object.data.materials.append(material)



def init_collections(collection_names=None, parent_collection=None):
    '''
    Make new collections or clear existing ones from objects (only touching
    collections specified in COLLECTION_NAMES, module level variable).
    
    Attributes
    ----------
    collection_names : list of strings
        If given, create or clear collections with the given names. Otherwise use
        module level variable COLLECTION_NAMES.
    parent_colelction : string
        If given, link to the collection with this name (instead Scene Collection)
    '''
    if collection_names == None:
        _collection_names = COLLECTION_NAMES
    else:
        _collection_names = collection_names

    for name in _collection_names:
        print("Clearing collection {}".format(name))
        collection = bpy.data.collections.get(name)
        
        if type(collection) != type(None):
            for obj in collection.objects:
                bpy.data.objects.remove(obj, do_unlink=True)
            
            # delete collection
            bpy.data.collections.remove(collection)
            
        # Then create new collection
        collection = bpy.data.collections.new(name)
        
        if parent_collection == None:
            bpy.context.scene.collection.children.link(collection)
        else:
            parent_collection.children.link(collection)
    
        # Shortcut dict for collections
        COLLECTIONS[name] = collection



def create_rhabdomeres(individual_objs=False, mirrors=True, R78=False, ommatidia_type=''):
    '''
    Create the object plane primitives according to the OBJECT_PLANE_TYPE constant.
    
    
    individual_objs : bool
        Create rhabdomeres as individual objects
    mirrors : bool
        Generate mirror pairs for left/right and frontal/dorsal
    R78 : bool
        If true, create separate R7 and R8. Othewsise, merge.
    ommatidia_type: string
        'pale' or 'yellow' or any aribtrary string. Only affects on the suffixes
        of the created primitives.
    '''
    if params.OBJECT_PLANE_TYPE == 'rhabdomeres':
        
        def create_rhabdomeres(mirror=False, mirror_lr=False):
            
            locations = params.RHABDOMERE_LOCATIONS.tolist()
            diameters = params.RHABDOMERE_DIAMETERS.tolist()
            if R78:
                locations.append(copy.copy(locations[-1]))
                locations[-2][2] = params.RHABDOMERE_DEPTH/4
                locations[-1][2] = -params.RHABDOMERE_DEPTH/4
                
                diameters.append(diameters[-1])
            
            for i_rha in range(len(locations)):
                loc = locations[i_rha]
                
                depth = params.RHABDOMERE_DEPTH
                loc[2] = loc[2] - depth/2
                
                if mirror:
                    loc = (loc[0], -loc[1], loc[2])
                if mirror_lr:
                    loc = (-loc[0], loc[1], loc[2])
                
                if R78 and i_rha >= 6:
                    depth /= 2
                    
                bpy.ops.mesh.primitive_cylinder_add(vertices=RHABDOMERE_VERTICES, radius=diameters[i_rha]/2,
                    depth=depth, location=loc)
                
                obj = bpy.context.active_object
                obj.name = "single_rhabdomere_{}".format(i_rha)
                bpy.ops.collection.objects_remove_all()
                COLLECTIONS['primitives'].objects.link(obj)
            
            
            if individual_objs:
                return None
            
            vertices_for_groups = {}
            
            def roundz(vector):
                return (np.round(1000000*vector)/1000000).tolist()
            
            # Join the rhabdomeres into a single object/mesh
            bpy.ops.object.select_all(action='DESELECT')
            for obj in COLLECTIONS['primitives'].objects:
                if "single_rhabdomere" in obj.name:
                    obj.select_set(True)
                    
                    vertices_for_groups[obj.name] = [roundz(v.co+obj.location) for v in obj.data.vertices]
                    
            bpy.context.view_layer.objects.active = obj
            bpy.ops.object.join()
            
            # Add vertex groups for each rhabdomere
            for group_name, vertices in vertices_for_groups.items():
                group = obj.vertex_groups.new(name=group_name)
                
                indices = [v.index for v in obj.data.vertices if roundz(v.co+obj.location) in vertices]
                
                group.add(indices, 1, 'ADD' )

            # Add materials
            materials = ['inactive', 'Rh1', 'Rh3', 'Rh4', 'Rh5', 'Rh6']
            for mat in materials:
                set_material(obj, mat, append=True)

            bpy.ops.transform.rotate(value=radians(-90), orient_axis='X', center_override=(0,0,0))
            
            # Set name for the joined rhabdomere
            objname = "rhabdomeres"
            if mirror:
                objname += '_mirror'
            if mirror_lr:
                objname += "_mirrorlr"
            if ommatidia_type:
                objname += '_{}'.format(ommatidia_type)
            obj.name = objname
            set_material(obj, "rhabdomere")
            
            if params.RHABDOMERE_ROTATION:
                bpy.ops.object.mode_set(mode='EDIT', toggle=False)
                bpy.ops.transform.rotate(value=params.RHABDOMERE_ROTATION, orient_axis='Y')
                bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
                
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        
        create_rhabdomeres(mirror=False)
        if mirrors:
            create_rhabdomeres(mirror=True)
            create_rhabdomeres(mirror=False, mirror_lr=True)
            create_rhabdomeres(mirror=True, mirror_lr=True)
    
    
        def create_hexagon():
            bpy.ops.mesh.primitive_cylinder_add(vertices=6, radius=1.0/2, depth=params.RHABDOMERE_DEPTH,
                location=(0,-params.RHABDOMERE_DEPTH-params.RHABDOMERE_DEPTH/2,0))
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.transform.rotate(value=params.HEXAGON_ROTATION, orient_axis='Z')
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            
            bpy.ops.transform.resize(value=params.HEXAGON_DIMENSIONS)
            bpy.ops.transform.rotate(value=radians(-90), orient_axis='X')
            
            obj = bpy.context.active_object
            obj.name = 'hexagon'
            bpy.ops.collection.objects_remove_all()
            COLLECTIONS['primitives'].objects.link(obj)
            set_material(obj, "Dark")
        
        create_hexagon()
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')



def create_lens(merge_surfaces=False, close_gap=False):
    '''
    Creates the facet lens primitive by calling the create_lens.py script
    with EXTERNAL_PYTHON.
    
    Adds objects named "lens_inner" and "lens_outer" into the "primitives" collections.
    
    OPTIONS
    ---------
    merge_surfaces  If true join the outer and inner lens surfaces together for better performance
                    In this case, only one object with name "lens" is added to the "primitives" collection
    close_gap       (experimental) If true, close the cap between the surfaces. Does not work with materials atm.
    '''
    # Call system python to run our lens making script (outside the bundled Blender Python)
    args = [EXTERNAL_PYTHON, os.path.join(MODELPATH, 'create_lens.py'),
        str(params.LENS_R), str(params.LENS_T), str(params.LENS_D), str(LENS_DXDY),
        str(close_gap), os.path.join(MODELPATH, 'tmp')]
    
    # Use previous lens data if parameters not changed
    args_cache = os.path.join(MODELPATH, 'tmp', 'lens_args.npy')
    try:
        with open(args_cache, 'r') as fp: past_args = fp.read()
    except:
        past_args = None
    
    if str(args) == past_args:
        print('Using lens mesh from past run (not changed)')
    else:
        subpr = subprocess.run(args)
        with open(args_cache, 'w') as fp: fp.write(str(args))
    
    # Load the npy file created by create_lens.py (to be run with Python external to
    # Blender; by default, Blender Python has no scipy package)
    vertices = np.load(os.path.join(MODELPATH, 'tmp', 'lens_vertices.npy'))
    faces = np.load(os.path.join(MODELPATH, 'tmp', 'lens_faces.npy'))
    
    if close_gap:
        iterover = zip([''], [90])
    else:
        iterover = zip(['_outer', '_inner'], [90, -90])
    
    for name, rotation in iterover:
        mesh = bpy.data.meshes.new("lens"+name)
        obj = bpy.data.objects.new(mesh.name, mesh)
        obj.name = "lens"+name
        col = bpy.data.collections.get("primitives")
        col.objects.link(obj)
        bpy.context.view_layer.objects.active = obj
        
        mesh.from_pydata(vertices.tolist(), [], faces.tolist())
        
        # Translate and rotate lens to right place
        bpy.ops.object.select_all(action='DESELECT')
        obj.select_set(True)
        bpy.ops.transform.translate(value=[0,params.LENS_TO_TIPS_DISTANCE,0])
        bpy.ops.transform.rotate(value=radians(rotation), orient_axis='X')

        # FIXME material setting like this doesnt't work when the close_gap is set
        set_material(obj, "Glass"+name)
        obj.select_set(True)
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
        
    if merge_surfaces and not close_gap:
        bpy.data.collections['primitives'].objects['lens_inner'].select_set(True)
        bpy.data.collections['primitives'].objects['lens_outer'].select_set(True)
        bpy.ops.object.join()
        bpy.context.view_layer.objects.active.name = "lens"
    
    
    bpy.ops.object.select_all(action='DESELECT')
    

    
def create_screening_pigments():
    '''
    Creates the screening pigments primitive that optically isolate
    the neighboring ommatidia from each other.
    '''
    cylinders = []
    for ssize in [1, params.HEXAGON_INOUT_RATIO]:
        bpy.ops.mesh.primitive_cylinder_add(vertices=params.HEXAGON_VERTICES,
            radius=ssize*1.0/2, depth=params.LENS_TO_TIPS_DISTANCE/ssize)
        bpy.ops.object.mode_set(mode='EDIT', toggle=False)
        bpy.ops.transform.rotate(value=params.HEXAGON_ROTATION, orient_axis='Z')
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        bpy.ops.transform.resize(value=params.HEXAGON_DIMENSIONS)
        bpy.ops.transform.rotate(value=radians(-90), orient_axis='X')
        bpy.ops.transform.translate(value=[0,params.LENS_TO_TIPS_DISTANCE/2-params.RHABDOMERE_DEPTH,0])
        
        obj = bpy.context.active_object
        obj.name = 'pigments'+str(ssize)
        bpy.ops.collection.objects_remove_all()
        COLLECTIONS['primitives'].objects.link(obj)
        
        cylinders.append(obj)
        
        bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
    

    
    bool = cylinders[0].modifiers.new(name='bool', type='BOOLEAN')
    bool.object = cylinders[1]
    bool.operation = 'DIFFERENCE'
    bpy.ops.object.modifier_apply(
            {"object": cylinders[0]},
            apply_as='DATA',
            modifier=bool.name)
    
    # Remove inner piece
    bpy.ops.object.select_all(action='DESELECT')
    cylinders[1].select_set(True)
    bpy.ops.object.delete()
    
    set_material(cylinders[0], "Dark")
    

def ommatidia(x,y,z, mirror=False, mirror_lr=False, add_to_collection='ommatidia',
        lens_normal=None, rhabdomere_movements=None, merge_stationary=False,
        ommatidia_types=[], lowlevel=True ):
    '''
    Create an ommatidia with lens, screening pigments, rhabdomere plane
    using the precreated objects (creates a linked copy).
    
    x,y,z               World x,y,z coordinates. +Y is the direction there lens normal point to
    mirrror             Hemisphere mirror rhabdomeres.
    mirror_lr           Left eye / right eye mirror rhabdomeres. If false create for right eye
    add_to_collection   Name of the collection where objects are added
    lens_normal         Lens normal vector to set the ommatidium orientation.
    rhabdomere_movements    In local rhabdomere coordinates, frame,x,y,z, a list of positions over time
    ommatidia_types : list of strings
        Randomly selects one of the specified ommatidia types
    lowlevel : bool
        If True, avoid bpy ops
    '''
    
    # Set nothing selected
    bpy.ops.object.select_all(action='DESELECT')
    
    rhabdomere_name = "rhabdomeres"
    if mirror:
        rhabdomere_name += "_mirror"
    if mirror_lr:
        rhabdomere_name += "_mirrorlr"
        
    if len(ommatidia_types) > 1:
        rhabdomere_name += "_{}".format(random.choice(ommatidia_types))
    
    if ONLY_RHABDOMERES:
        objnames = [rhabdomere_name, 'hexagon']
    else:
        if not 'lens' in bpy.data.collections['primitives'].objects:
            objnames = [rhabdomere_name, 'hexagon', 'lens_inner', 'lens_outer', 'pigments1']
        else:
            objnames = [rhabdomere_name, 'hexagon', 'lens', 'pigments1']
    
    #if merge_stationary:
    #    for objname in objnames:
    #        obj = bpy.data.collections['primitives'].objects[objname]
    #        obj.select_set(True)
    #    bpy.context.view_layer.objects.active = obj
    #    bpy.ops.object.join()
    #    bpy.context.view_layer.objects.active.name = 'merged_stationary'
    
    objs = []
    
    if lowlevel:
        for objname in objnames:
            obj = bpy.data.collections['primitives'].objects[objname].copy()
            #mesh_offset = obj.location.copy()
            obj.location += mathutils.Vector((x,y,z))

            objs.append(obj)
            # No linking, remember to do it outside this function if lowlevel
        
        #for obj in objs:
        #    bpy.data.collections[add_to_collection].objects.link(obj)
    else:
        for objname in objnames:
            obj = bpy.data.collections['primitives'].objects[objname].copy()
            #obj.name = objname + "_{},{},{}".format(x,y,z)
            bpy.data.collections[add_to_collection].objects.link(obj)
            
            
            obj.select_set(True)
            bpy.ops.transform.translate(value=[x,y,z])
            obj.select_set(False)
            objs.append(obj)
            

        if rhabdomere_movements is not None:
            # Animate rhabdomere object only
            obj = objs[0]
            obj.select_set(True)
            baseloc = copy.deepcopy(obj.location)
            
            for frame,a,b,c in rhabdomere_movements:
                bpy.context.scene.frame_set(int(frame))
                obj.location = copy.deepcopy(baseloc)
                bpy.ops.transform.translate(value=[a,b,c], orient_type='LOCAL')
                obj.keyframe_insert(data_path="location", index=-1)
        
        
        for obj in objs:
            obj.select_set(True)
        
        if lens_normal is not None:
            a,b,c = lens_normal
            xrot = -atan2(c,b)
            bpy.ops.transform.rotate(value=xrot, orient_axis='X', orient_type='GLOBAL')
            
            zrot = -atan2(-a,b)
            bpy.ops.transform.rotate(value=zrot, orient_axis='Z', orient_type='GLOBAL')
    
    
    return objs




def get_ommatidia_locations(side):
    '''
    Getting the locations of ommatidia on the eye
    
    side        "left" or "right" ("left" mirrors x coordinates)
    '''
    
    def curved_eye_coordinates():
        '''
        Get coordinates from an ico-sphere
        '''
        verts = []
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=params.EYE_RADIUS)
        obj = bpy.context.active_object
        for i in range(10000):
            try:
                verts.append(obj.data.vertices[i].co)
            except:
                break
        
        return np.array(verts)
            
    
    if params.CURVED_EYE_METHOD == 'icosphere':
        locations = curved_eye_coordinates()
    elif params.CURVED_EYE_METHOD == 'external':
        
        
        # Call system python to run our lens making script (outside the bundled Blender Python)
        args = [EXTERNAL_PYTHON, os.path.join(MODELPATH, 'curved_eye_coordinates.py'),
            str(params.EYE_RADIUS), str(params.EYE_HORIZONTAL_R), str(params.LENS_D/2),
            os.path.join(MODELPATH, 'tmp')]
        
        # Use previous lens data if parameters not changed
        args_cache = os.path.join(MODELPATH, 'tmp', 'eye_coordinates_args.txt')
        try:
            with open(args_cache, 'r') as fp: past_args = fp.read()
        except:
            past_args = None
        
        if str(args) == past_args:
            print('Using ommatidia locations from a past run (not changed)')
        else:
            subpr = subprocess.run(args)
            with open(args_cache, 'w') as fp: fp.write(str(args))
        
        locations = np.load(os.path.join(MODELPATH, 'tmp', 'curved_eye_coordinates.npy'))
    
    
    # No  and allow only X>0
    if params.EYE_LOWER_ANGLE:
        k = -1/tan(params.EYE_LOWER_ANGLE)
        #lowest_z = - params.EYE_RADIUS * cos(params.EYE_LOWER_ANGLE)
        locations = [p for p in locations.tolist() if k*p[0]<p[2] and p[0]>-params.EYE_EXTEND_INVARDS]
        locations = np.array(locations)


    # Limit the maximal number of ommatidia
    if len(locations) > params.OMMATIDIA_LIMIT:
        locations = locations.tolist()
        locations.sort(key=lambda x: x[0], reverse=True)
        locations = np.array(locations[0:params.OMMATIDIA_LIMIT])
    
    # Center on Y axis, for some reason isn't centered always?
    avglocy = np.mean([y for x,y,z in locations])
    locations = [[x,y-avglocy,z] for x,y,z in locations]
    
    if side == 'left':
        locations = [[-x,y,z] for x,y,z in locations]
        
    return locations




def build_eye(side, locations, movement_data=None, move_along_x=True, scale_factor=None, ommatidia_types=[],
        lowlevel=False):
    '''
    Main eye building function.
    '''
    ommatidia_objs = []
    
    name = side + 'eye'
    
    #if scale_factor is not None:
    #    locations = np.array(locations)
    #    locations = locations * scale_factor

    # Make sure all the objects are unselected at the start; Our function ommatidia()
    # requires this.
    bpy.ops.object.select_all(action='DESELECT')
    
    
    for i_loc, location in enumerate(locations):
        print("{}, ommatidia {}/{}".format(name, i_loc+1, len(locations)))

        mirror = params.MIRROR_HEMISPHERES and location[1]<0
        if side == 'left':
            mirror_lr = True
        else:
            mirror_lr = False

        objs = ommatidia(*np.array(location), mirror=mirror, mirror_lr=mirror_lr,
            add_to_collection=name, ommatidia_types=ommatidia_types, lowlevel=lowlevel)
        
        x,y,z = location / np.linalg.norm(location)#params.EYE_RADIUS
        
        try:
            phi = asin(z)
        except:
            phi = pi/2
        
        if lowlevel:
            for obj in objs:
                obj.rotation_euler[0] += phi
        else:
            bpy.ops.transform.rotate(value=-phi, orient_axis='X', orient_type='LOCAL')
        
        try:
            phi = atan(x/y)
        except:
            phi = pi/2
        if y < 0:
            phi = phi + pi
        
        if lowlevel:
            for obj in objs:
                obj.rotation_euler[2] += -phi
                
                if move_along_x:
                    if side == 'right':
                        obj.location[0] += params.EYE_HEAD_CP_DISTANCE
                    else:
                        obj.location[0] -= params.EYE_HEAD_CP_DISTANCE
            
        else:
            bpy.ops.transform.rotate(value=phi, orient_axis='Z', orient_type='GLOBAL')
            
            
            # Rotation of the rhabdomeres across the eye
            [obj.select_set(False) for obj in objs]
            objs[0].select_set(True)
            bpy.context.view_layer.objects.active = objs[0]
            bpy.ops.transform.rotate(value=pi/2-atan2(y,x), orient_axis='Z', orient_type='LOCAL')
            bpy.context.view_layer.objects.active = None
            objs[0].select_set(False)
            
        
            if move_along_x:
                # Finally, move the eye in the X-axis to the right distance from the head center point
                [obj.select_set(True) for obj in objs]
                if side == 'right':
                    bpy.ops.transform.translate(value=(params.EYE_HEAD_CP_DISTANCE, 0.0, 0.0))
                else:
                    bpy.ops.transform.translate(value=(-params.EYE_HEAD_CP_DISTANCE, 0.0, 0.0))
                [obj.select_set(False) for obj in objs]


        
        if movement_data is not None:
            obj = objs[0]
            exaggerate = 10*scale_factor
            vector = movement_data[i_loc]
            spoint = np.array(obj.location)
            scene = bpy.context.scene
            
            for i_repeat in range(0,40, 1):
                scene.frame_set(i_repeat*20)
                obj.location = spoint
                obj.keyframe_insert(data_path="location", index=-1)
                
                scene.frame_set((i_repeat*20)+4)
                newloc = (spoint[0]+exaggerate*vector[0], spoint[1]+exaggerate*vector[1],spoint[2]+exaggerate*vector[2])
                obj.location = newloc
                obj.keyframe_insert(data_path="location", index=-1)
        
        ommatidia_objs.append(objs)
        

        
    return name, ommatidia_objs



def get_movement_data(side, locations):
    '''
    Side "left" or "right"
    '''
    
    # Save the locations in a npy file in the tmp folder
    locations_savefn = os.path.join(MODELPATH, 'tmp', 'locations.npy')
    np.save(locations_savefn, np.array(locations, dtype=np.float).tolist(), allow_pickle=False)

    subpr = subprocess.run([EXTERNAL_PYTHON, os.path.join(MODELPATH, 'animation_data.py'),
        side, locations_savefn, MODELPATH])
    
    animations = np.load(os.path.join(MODELPATH, 'tmp', 'movement_data.npy'))
    
    return animations

    
    
def hide_primitives():
    for obj in bpy.data.collections['primitives'].objects:
        obj.hide_render = True
        obj.hide_viewport = True
        

MODELPATH = bpy.path.abspath("//")
COLLECTION_NAMES = ['primitives', 'rhabdomeres', 'fakeeye', 'ommatidia', 'lefteye', 'righteye', 'head', 'jounisim']
COLLECTIONS = {}
params = getattr(__import__(__name__), EYE_DESIGN)



def main(only_ommatidium=False, pale_yellow=False, merge_lens_surfaces=False,
        animate=False):
    '''
    Initializes the eye by
        - deleting previous any previous ()
    
    only_ommatidium : bool
        If True, create only the ommatidium elements (in primitive collection)
    pale_yellow : bool
        If True, create separate objects for pale and yellow ommatidia,
        and distribute them randomly. Makes it easy to color them separately after.
    merge_lens_surfaces : bool
        If True, close the gap between the lens inner and outer surfaces (not perfectly atm.)
    '''
    start_time = time.time()
        
    bpy.ops.object.select_all(action='DESELECT')
    rescale1(SCALE_FACTOR)
    
    init_collections()
    bpy.ops.outliner.orphans_purge()
    
    if pale_yellow:
        ommatidia_types = ['pale', 'yellow']
    else:
        ommatidia_types = ['']
    
    for ommatidia_type in ommatidia_types:
        create_rhabdomeres(individual_objs=False, mirrors=True, R78=True, ommatidia_type=ommatidia_type)

    create_lens(merge_surfaces=merge_lens_surfaces)
    create_screening_pigments()
    
    if only_ommatidium:
        return None
    
    all_ommatidia = {}
    for side in ['left', 'right']:
        locations = get_ommatidia_locations(side)
        
        if animate:
            movement_data = get_movement_data(side, locations)
        else:
            movement_data = None
            
        collection_name, ommatidia = build_eye(side, locations, movement_data=movement_data,
            move_along_x=True, scale_factor=SCALE_FACTOR, ommatidia_types=[''],
            lowlevel=LOWLEVEL_OPS)
        
        all_ommatidia[collection_name] = ommatidia


    if LOWLEVEL_OPS:
        # Lowlevel: Linking everything to the scene as late as possible
        for side in all_ommatidia:
            for objs in all_ommatidia[side]:
                for obj in objs:
                    bpy.data.collections[side].objects.link(obj)
    
    
    hide_primitives()
    
    print("FINISHED!")
    print('Total initialization time: {} minutes'.format((time.time()-start_time)/60.))
    
    
if __name__ == "__main__":
    main()
