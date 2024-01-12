# Block-Modeling-Matlab-Codes
Block modeling codes to solve for fault slip rates and block rotations from geodetic data

These codes have been used by scientists and students at the University of Nevada, Reno to analyze crustal deformation in tectonic plate boundary zones with complex faulting. 

The codes are provided as scripts executable in Matlab.  It is necessary to have access to the Matlab platform in order to use them. The Mapping toolbox code distance.m is used in some places but can be replaces with the baz.m code (also provided) if you do not have access to the Mapping toolbox. 

To install the codes simply place the files in a folder.  Make sure your Matlab paths can see those folders (some "addpath" commands are included in the example script to help with this).  

An example script is provided that illustrates a typical workflow.  It reads in a model and a GPS velocity field, makes plots, solves the block rotation/slip rate equations, and outputs results in tables and graphics. If all is working OK the example script should run to completion without errors. See comments within the codes for some detailed information.

That manuscript describing analysis of crustal deformation in the western Great Basin using these codes can be cites as:
Hammond, W. C., G. Blewitt, and C. Kreemer (2011), Block modeling of crustal deformation of the northern Walker Lane and Basin and Range from GPS velocities, J. Geophys. Res., 116, B04402, doi:10.1029/2010JB007817.  

-Bill Hammond,
Nevada Geodetic Laboratory,
University of Nevada, Reno
