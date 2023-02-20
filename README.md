# CG compound eye

*GG-compound-eye* is a parametric computer graphics (CG) model of
the *Drosophila* compound eyes, capable of producing the deep pseudopupil (DPP).

It uses Blender's Python-scripting capabilities to procedurally build the model,
and LuxCoreRender to simulate/render.


## Prerequisites

To produce the eye 3D model

* Blender 2.80, or newer (tested with 2.93)
* Working Python 3 installation with Scipy >=1.6.0

Additionally, to render the deep pseudopupil

* Luxcorerender >= 2.3 (tested with 2.3)


## How to use

1) Ensure the files here and your Blender File (.blend) are in the same directory
1) In Blender's Text Editor (Scripting), open *model_init.py*
1) Press *Run Script* button

This leads to the 3D model initialization taking some minutes.
The model is quite small so you have to zoom in (it fits inside the default cube).

Process and possible error messages can be monitored from the system console
(On Windows, *Window -> Toggle System Console*; On other platforms, launch Blender from a shell)

To render DPP images
* configure materials (not done by the script yet)
* change the rendering engine to LuxCore
* enable camera's Depth of Field


## Modification

Distributed as a script, the CG-eye can be easily modified.
For details, see the documentation in *model_init.py*.


## Performance

The model building blocks (lenses, rhabdomeres, ...) are created in the beginning
(see the hidden collection *primitives*), and in each ommatidia,
the corresponding objects merely reference this data.
These *linked duplicates* keep the memory usage relatively low (~220 MiB).

However, since each ommatidia has 5 parts, already with 700 ommatidia there are
3500 objects in the scene. This amount of objects already slows
Blender in some operations.
For example, deleting the model can hang up for some time.

I have optimized the model creation using the Blender's low-level API and linking the objects
to the scene as late as possible (see LOWLEVEL_OPS, turned on by default).
This strategy allows creating larger models up to some thousand ommatidia
relatively fast.
