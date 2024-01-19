Readme file for codes that build a block model from fault data

Generating a block model from raw fault inputs can save a lot of time in building models which would otherwise be built 'by hand' from fault data.  These codes were first developed as a way to reduce labor in generating block models and were then later incorportated into iterative procedures that build thousands of models to improve statistics of slip rate estimation in zones of complex fault networks.  The core framework and reference for the block modeling approach is given in the other README.md file for this GitHub project. This file "README_buildblockmodelcodes.md" is for just the part of the codes that takes fault data and builds a block model, i.e. generates Matlab variables containing lists of node coordinates 'nodes', block boundaries specifications as ordered lists of nodes 'blocks', and faults segments that have locking depths and dips 'faults'.  

The reference for these codes is
Hammond, W.C., C. Kreemer, and G. Blewitt, 2024, Robust imaging of fault slip rates in the Walker Lane and Western Great Basin from GPS data using a multi-block approach, Journal of Geophysical Research-Solid Earth, in revision. 

To learn about and test the codes use the example "ExampleMakeBlocks.m" that is inside the buildblockmodelcodes folder.

The procedure has several steps:

1) Importing the fault data.
Fault data must be brought into a Matlab structured array variable "geol".  This can be the most time consuming part of the process, and these codes do not perfectly automate it.  However, an example is provided in the script "ReadUSGSFaults.m" to show how it can be done.  In this case data from the USGS National Seismic Hazard Modeling Project database was used, obtained from shapefiles provided by the USGS. The script contains the URL from which data were obtained. The Matlab script reads in geologic fault data from the USGS shapefile, including latitude/longitude for points along the trace, dip, locking depth, rake and other information about the faults.  That cose develops a set of faults that is truncated to include only those within the area of interest.  

Once the USGS fault data is read in, there is then some work checking the faults to find any that have segments that cross other segments. Crossing segments are not allowed in this modeling framework so those must be identified and removed.  Read the code "ReadUSGSFaults.m" to see how crossing segments are identified and removed.  You many have to run the code multiple times to identify bad faults and specify them for removal. In the end a Matlab workspace file is saved that has the variable "geol" and we are then ready for the next step.

2) Build a block model
The codes in the folder "buildblockmodelcodes" are used to generate a set of nodes and block boundaries that represent the fault network.  The algorithm starts with the data on fault traces, simplifies them with approximately evenly spaced segments, seeds the boundary with nodes and then forms a Delaunay triangulation of them.  The resulting tesselation by triangles can be thought of as a primitive block model with a large number of small blocks. 

The strategey is then to systematically reduce the number of blocks by combining adjacent blocks into larger blocks, while at each step making sure the block model is still well posed, maintains block contiguity and boundaries along the fault segments, and carries forward the dip/locking depth information about the faults. A set of codes hunts for candidate block combinations based on user adjustable parameters that specify maximum block number, minimum block size, aspect ratio, convexity.  A scoring system for block bounding segments helps place controls on which can be eliminated without making an  block model invalid. There are checks along the way that make sure model integrity is maintained.   

The algorithm iterates until one of several criteria are reached. Usually at that point no further improvements can be made given the selected paramters. Iterations cease and the variables 'nodes','blocks','faults' are save to a Matlab workspace for the next step.  If the verbose option is selected many plots will show the progress of block model reduction as the iterations proceed.   

3) Solve for slip rates
The process in most cases generate a 'legal' block model which can then be read into a script that runs the code 'BlockInter.m' which solves for slip rates and block rotations.  The codes to do that are in the other folder "blockcodes". In the folder "Example" there is a script "RunBlock.m" which runs through the workflow of solving for slip rates and block rotations.  If you want to run the model you built here, make sure you load it in by specifying the right path and file name in "RunBlock.m".

Try to run the example scripts "ExampleMakeBlocks.m" and "RunBlock.m" first.  If you can get those to work with no errors then you should be able to change the area of interest and use different fault data to construct your model.  


