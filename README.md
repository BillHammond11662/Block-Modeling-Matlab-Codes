Readme file for Block-Modeling-Matlab-Codes

Fault slip rates are key parameters that feed into seismic hazard analyses. An increasingly accepted method for estimating fault slip rates is to use a geodetic block model, where fault data and velocities of active crustal motion are combined analytically to solve for the relative motion between crustal blocks over long periods of (geologic) time.

However, some active plate boundary zones have many faults that form complex systems that make model construction difficult and non-unqiue.  The codes presented here form a tool kit for building, solving, and plotting results of block models of active crustal deforamtion. The approach here is to provide modular software, giving the user flexibility in how to construct and change the models, and regularize the solutions.  These codes have been used by a number of academic scientists and students to analyze active crustal deformation in tectonic plate boundary zones with complex faulting. The codes were originally developed to solve for slip rates and crustal block rotations in the Northern Walker Lane, western Great Basin of the United States.  

The codes are provided as scripts executable in Matlab.  It is necessary to have access to the Matlab platform in order to use them. The Mapping toolbox code "distance.m" is used in some places but can be replaced with the code baz.m (also provided) if you do not have access to the Mapping toolbox. Otherwise no specialized Matlab toolboxes are needed. 

To install the codes simply place the files in a folder on your computer.  Make sure your Matlab paths can see those folders (some "addpath" commands are included in the example script to help with this, though you may have to adjust them depending on how you store folders).  

An example script "RunBlock.m" in the folder "Example" is provided that illustrates a typical workflow.  It reads in a model and a GPS velocity field, makes plots, solves the block rotation/slip rate equations, and outputs results in tables and graphics. If all is working OK the example script should run to completion without errors. See comments within the codes for some detailed information.

That manuscript describing the analysis of crustal deformation in the western Great Basin using these codes is:
Hammond, W. C., G. Blewitt, and C. Kreemer (2011), Block modeling of crustal deformation of the northern Walker Lane and Basin and Range from GPS velocities, J. Geophys. Res., 116, B04402, doi:10.1029/2010JB007817.  
which can be cited if you are publishing using these codes.

Also see the following for examples of their use:

Hammond, W.C., C. Kreemer, and G. Blewitt, 2024, Robust imaging of fault slip rates in the Walker Lane and Western Great Basin from GPS data using a multi-block approach, Journal of Geophysical Research-Solid Earth, in revision. 

Young, Z, C. Kreemer, W.C. Hammond, and G. Blewitt, 2023, Interseismic strain accumulation between the Colorado Plateau and the Eastern California Shear Zone: Implications for the seismic hazard near Las Vegas, Nevada, Bulletin of the Seismological Society of America, 113(2), 856-876, https://doi.org/10.1785/0120220136.

Li, X, W.C. Hammond, I.K.D. Pierce, J. M. Bormann, Z. Zhang, C. Li, W. Zheng, P. Zhang, 2023, Present-day strike-slip faulting and intracontinental deformation of North China: Constraints from improved GPS observations, Geochemistry, Geophysics, Geosystems, v. 24(7), GGGE23102, https://doi.org/10.1029/2022GC010781.

Li, X., I.K.D. Pierce, J. Bormann, W.C. Hammond, C. Li, W. Zheng, P. Zhang, 2021, Contemporary deformation of the northeastern Tibetan Plateau and its surroundings revealed from GPS block model, J. Geophys. Res. - Solid Earth, v 126, 5, e2020JB020733, special section on the 100-year anniversary of the great 1920 Haiyuan Earthquake: What have we learnt on large continental earthquakes and faults?. https://doi.org/10.1029/2020JB020733.

Bormann, J., W.C. Hammond, C. Kreemer, G. Blewitt, 2016, Accommodation of missing shear strain in the Central Walker Lane, western North America: Constraints from dense GPS measurements, Earth and Planetary Science Letters, 440, 169-177, doi:10.1016/j.epsl.2016.01.015.


-Bill Hammond,

Nevada Geodetic Laboratory,
University of Nevada, Reno


