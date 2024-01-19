% Example use of MakeBlockModel.m
%
% Generates a block model from fault data.
% Makes "blockmodel.mat" which can be used as input to block modeling codes
% provided on this site.
%
% Please see README_buildblockmodelcodes.txt for a brief description of
% these codes.
% 

clear
close all;  % Good to delete all plots to start.  Otherwise desktop gets crowded. 

tic;        % starts the run time clock

% The diary records text output to the Matlab command window.
% This is useful for seeing steps of model generation
delete('MakeBlockModelExample.txt');
diary MakeBlockModelExample.txt;
diary on;

% make sure the folders with codes are in your Matlab paths.  These are
% examples, not necessarily how you will organize on your system...
addpath('../blockcodes');
addpath('../blockcodes/figurecodes');


%% Get the input fault data from the database
%
% This section loads the data structure that contains the fault data: 'geol'
% these data were obtained from the USGS National Seismic Hazard Mapping
% Database.  A section was clipped out for the purpose of this example. 
%
% Hatem, A.E., et al., 2021, Earthquake geology inputs for the National Seismic 
% Hazard Model (NSHM) 2023, version 1.0: U.S. Geological Survey data release, https://doi.org/10.5066/P918XCUU.
% obtained at https://www.sciencebase.gov/catalog/item/5fe3bbbcd34ea5387deb4a85
%
% read the script ReadUSGSFaults.m to see how we made the variable 'geol'
% from the USGS fault database. 

load('GeolForGitHubExample.mat','geol');
Nf = length(geol);

% Set the lat/lon boundary of the model
lonmin=-119;
lonmax=-117;
latmin=35.5;
latmax=37.5;
bounds = [lonmin lonmax latmin latmax];

% Make a matlab plot of the input fault network.
% The box shows area that will be modeled with blocks

figure(1);
clf;
PlotStates(0,0);
hold on;
plot([lonmin lonmax lonmax lonmin lonmin],[latmin latmin latmax latmax latmin],'k-');
axis(bounds+[-.2 .2 -.2 .2]);
set(gca,'DataAspectRatio',[1 cosd(mean([latmin;latmax])) 1]);
for i=1:Nf
    plot(geol(i).lonseg,geol(i).latseg,'r-','linewidth',2);
end
drawnow;

%% Perform the iterative block model construction

% First set some flags that control behavior:

optverb = 1;  % If optverb==1 more text comments are provided in Matlab command window

optplot = 1;  % Warning: If optplot==1 this generates a lot of matlab plots 
              % that provide detail about the iterative construction of the model.
              % It is ok to set this to 0 after seeing how the reduction of
              % block numbers by sequential combination of smaller blocks into larger blocks works. 
              % Seeing the iterations in detail can help with debugging.

optdamp = 1;    % Note: If optdamp==0 then the only block boundaries that are given fault elements 
                % are those that are bounded by input fault geometries present in the input fault
                % database 'geol'.  If optdamp==1 then ALL block boundaries,
                % regardless of whether an input fault is present in 'geol'
                % When fault elements are present.
                % For all segments in the variable 'faults' model parameters 
                % for fault slip rate are generated in BlockInter.m, which allows for damping of the solution
                % for slip rate at these locations.  Thus when optdamp==1
                % it is possible to regularize relative motion between blocks
                % ('slip') even if there is no real fault there in the
                % database.  This is useful for model regularization and
                % perhaps damping on-fault vs. off-fault deformation
                % separately.

disp(' ');
disp('Building Block Model Geometries');


% MakeBlockModel.m creates a block model from the fault data.
% It starts by making a model with a lot of triangular blocks from a
% Delaulay triangulation of the nodes that trace the faults and that lie on the bounding lat/lon box.
% The algorithm then systematically reduces the number of
% blocks to make the model simpler and better conditioned before inverting
% for slip rates.  At each step the algorithm attempts to combine blocks
% together in an intelligent way.  Blocks that are small, are very narrow, have concavities, 
% are favored for removal.  Block boundaries that lie on faults or are along 
% the edge of the bounding lat/lon box are preserved.
%
% The following parameters can be adjusted to help guide the iterative
% refinement of the model as the algorithm works to combine blocks to
% reduce their number.  

Bnmax = 50;         % maximum number of blocks, will quit after getting below this (but may quit soon if it runs out of options).
IntAngleMax = 190;  % maximum allowed interior block polygon angle. If >180 indicates concavity in block.
IntAngleMin = 20;   % minimum allowed interior block polygon angle
Amin = 170;         % minimum block area in km^2
Rmax = 1.9;         % maximum shape parameter allowed (circle=1)
Degsmall=1;         % minimum interior angle for sliver blocks on boundary
params={Bnmax;IntAngleMax;IntAngleMin;Amin;Rmax;Degsmall};

% Most of what happens goes on in MakeBlockModel.m
[nodes,blocks,faults,geoltrunc,ifltdamp]=...
    MakeBlockModel(geol,bounds,params,optverb,optplot,optdamp);

%% Print out some facts about the model geometries and make some plots

% blockgeoms.m provides geometric charcteristics of the model blocks
[A,R,T,C,P,amax,amin,Nv] = blockgeoms(blocks,nodes);

Nb = length(A);
disp(' ');
disp('Block   N       N      Area  Perim.   Shape    InteriorAngles ');
disp('  #   Verts  Concave  (km^2)  (km)    Param.   Min(˚)  Max(˚) ')

for i=1:Nb  
    str = [sprintf('%3.0f',i) '    ' sprintf('%2.0f',Nv(i)) '      ' sprintf('%2.0f',C(i))  '     ' sprintf('%6.1f',A(i))  ' ' sprintf('%6.1f',P(i)) ' ' sprintf('%6.1f',R(i)) ...
       '    ' sprintf('%6.1f',amin(i)) ' ' sprintf('%6.1f',amax(i))  ];
    disp(str);
end

PlotFinalBlockModel;

save('blockmodel','nodes','faults','blocks');

diary off;

toc

