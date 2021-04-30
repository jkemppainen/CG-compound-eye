# CG compound eye

*GG-compound-eye* is a parametric computer graphics (CG) model of
the *Drosophila* compound eyes, capable of producing the deep pseudopupil (DPP).

It uses in-Blender-scripting to initialize the 3D model
and LuxCoreRender to simulate/render.


## Prerequisites

* Blender >= 2.8 (tested with 2.8.1)
* Luxcorerender >= 2.3 (tested with 2.3)
* Working Python 3 installation with Scipy >=1.6.0

## How to use

1) Ensure the files here and your Blender File (.blend) are in the same directory
1) In Blender's Text Editor (Scripting), open *model_init.py*
1) Press *Run Script* button

This leads to the 3D model initialization taking some minutes.
Process and possible error messages can be monitored from the system console
(On Windows, *Window -> Toggle System Console*; On other platforms, launch Blender from a shell)

To render DPP images
* configure materials (not done by the script yet)
* change the rendering engine to LuxCore
* enable camera's Depth of Field

## Modification

Distributed as a script, the CG-eye is easily modified.
For details, see the documentation in *model_init.py*.

## Troubeshooting

Make sure that

* `3D cursor` is at the world origin
* Blender is in the `Object Mode`
* `EXTERNAL_PYTHON`'s value, when typed on a shell (default *python*), opens a Python interpreter where `import scipy` works

Try also with a new blend file in case of problems.
