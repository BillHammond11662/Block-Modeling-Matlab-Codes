% Run a Block Model
%
% This loads into Matlab the variables 'nodes', 'blocks', and 'faults' that represent 
% the geometry of the blocks, including the coordinates of the 'nodes' that define the 
% perimeters of 'blocks' which contiguously fill the space, and 'faults' which lie along some block
% boundaries.  Not all block boundaries must have faults.  Only block boundaries that have faults
% will have elastic strain accumulation owing to deep dislocations incorporated into the model
% and have slip rates solved for. 

clear;
close all;
tic

addpath('../blockcodes');
addpath('../blockcodes/figurecodes');


%%  Get model geometries. Blocks and faults. Make some plots.

% The matlab workspace file 'blockmodel.mat' has variables 'nodes', 'blocks',
% and 'faults'.  These are the necessary variables to define the block model
% geometries. 
%
% 'nodes' is Nx2. Where N is number of nodes in the entire model. 
%                 nodes(:,1) are longitudes, nodes(:,2) are latitudes 
% 'blocks' is a structured array.  
%                 blocks.block1 is a list of nodes that defines the perimeter of block1.
%                 There can be any number of blocks, but ideally there are
%                 no holes in blocks and where they abut one other the
%                 nodes of the boundaries are listed in both the blocks
%                 that touch one another.
% 'faults' is Px4 where P is the number of faults. 
%                   faults(:,1) and faults(:,2) are the nodes numbers that define
%                   the line segment which is the surface trace of the
%                   fault.  
%                   faults(:,3) are the locking depths of each fault in km
%                   faults(:,4) are the dips of the faults in degrees.


disp('Reading in block model from the .mat file.')
load blockmodel.mat;
% If you want to run the model generated by codes in the folder
% "buildblockmodelcodes" comment out the previous line and uncomment the
% next line:
% load('../buildblockmodelcodes/blockmodel.mat');
% Note the sample code below in the file SimplifyBlocks.m only works on
% the version of 'blockmodel.mat' that is in the folder "Example".

% Uncomment/comment out the next line to run or not run the code SimplifyBlocks.m.
% This codes only works on the specific version of blockmodel.mat that is
% in the folder "Example". SimplifyBlocks.m is provided as an example of
% how blocks can be refined 'by hand' by choosing block boundaries to
% eliminate.  You can use PlotBlocks.m (with opt3==1) to see node numbers of segments that
% you wish to eliminate.
%
% SimplifyBlocks;  % uncomment this line to 'simpify' the model.

disp('Performing Model Checks...')
CheckNodes(nodes, faults, blocks);
% Note that "CheckNodes" gives a warning that node 98 is very close to node 99.  
% This is not fatal, just a warning, so is OK in this case to forge ahead

CheckFaults(nodes,blocks,faults);

bounds = [min(nodes(:,1)) max(nodes(:,1)) min(nodes(:,2)) max(nodes(:,2))] + [-.2 .2 -.2 .2];

Figure1;
Figure2;
Figure3;


%%  Get GPS velocities that constrain block motion.  
% These are from the MIDAS velocity field in North America reference frame
% available at http://geodesy.unr.edu/velocities/midas.NA.txt
% a couple of outlier velocities were removed

velfile = 'example_midas_vel.txt';
[lat,lon,vn,ve,vu,sn,se,su,corrne,correu,corrnu,sta]=ReadVelFile(velfile);

% Plot them.
Figure4;
Figure5;


%% Run the block model, i.e., solve for block rotation and slip rates
% using blocks,nodes,faults to definte block geometery and 
% ve,vn are velocities with 
% se,sn are uncertainties in velocity
% vu,su are not used.

blocknames = fieldnames(blocks);
M=size(blocknames,1);
P=size(faults,1);

% specify some basic parameters for the block model
nu = .25;           % Poisson's ration of crust (unitless)      
mu = 3e10;          % shear modulus of crust (Pa)

% These regularize the model by specifying weights or damping for each model parameter
% Regularization provides analytical control over the solution, e.g., what
% should happen if there is not enough data to constrain block motion,
% should the spin be minimized or the sum of slip rates on the block
% perimeter?  Sometimes there are not enough GPS stations on a block so some
% assumption must be made. 
%
% Play with these to see what effects they have on the slip rates and
% block rotation rates.
alpha = 0.00001;    % apriori uncertainty in slip rate consistency constraint, units in "m/yr"
                    % alpha should be small to enforce kinematic
                    % consistency between block rotations and slip rates.
                    % I don't recommend changing this value.
beta = 3e-9;        % apriori uncertainty in omega in units "rad/yr" for weighting
                    % small values allow for less rotation.
gamma = .0020;      % apriori uncertainty in slip rates in units "m/yr" for weighting
                    % small values allow for less fault slip
vabeta = 1;         % if vabeta==1 then vertical axis spin is damped instead of total rotation Euler pole
                    % this is an important distinction because blocks may
                    % translate rapidly (have an Euler pole 90˚ distance
                    % from block) and have very little vertical axis spin (rotation
                    % around axis pointing towards center of block).  So
                    % with vabeta=1 you damp the spin rate, not the translation
                    % rate.
L=6;                % L is max number of nearest fault's interseismic strain from locking 
                    % contributing to each station's predicted velocity. 

% If creep=[] we will not account for active creep on faults
creep=[];
creep_a=[];

%
% slip_a is the apriori value of slip rate for each fault. Nans where none
% is assumed.
% s_slip_a is the uncertainty on the assumed slip_a.  
% If slip_a=[] apriori values of slip rates are not applied for any fault.  
% However, if there are constraints from geologic fault slip rates you can use
% these as initial constraints on the model.  For example if fault #137 is
% thought to have a strike slip rate of 2.0±0.2 mm/yr dextral and a normal slip rate of 0.0±0.2 mm/yr you could put: 
% slip_a(131:137,1)=-.002;
% s_slip_a(131:137,1)=.0002;
% slip_a(131:137,2)=0;
% s_slip_a(131:137,2)=.0002;
%
% Note: conventions about slip rate signs are giving in help for BlockInter.m
% slip(:,1) are strike slip rates, slip(:,2) are dip slip rates (in along dip direction not projected to horizontal)
%
% But for now we don't assume any:
slip_a=nan(P,2);
s_slip_a=nan(P,2);

% BlockInter.m can also solve for uniform horizontal tensor strain rate
% within any block.  This is train on top of any coming from interseismic
% strain accumulation on blocks. This could be useful of you have a large
% block that is deforming uniformly.
%
% But for now apriori values of strain rates are not applied (so we put
% nans)
zeta = 0;
strain_a = nan(M,3);
s_strain_a = nan(M,3);

% Apriori values can also be imposed for block rotations if desired.
% But in this case they are not applied (so we put
% nans)
omega_a = nan(M,3);
s_omega_a = nan(M,3);

% BlockInter.m builds a matrix to do a linear inversion for block
% rotation rates and fault slip rates simultaneously.
% The matrix is inverted and results passed to output variables.
tic
disp('  ')
disp('Solving for block rotations and strain rates')
[omega,s_omega,slip,s_slip,strain,s_strain,covm,vnpred,vepred,chi2,dof,ilist,mse]=...
    BlockInter(lon,lat,ve/1000,vn/1000,vu,se/1000,sn/1000,su,nodes,blocks,...
    faults,nu,mu,alpha,beta,gamma,zeta,L,omega_a,s_omega_a,slip_a,s_slip_a,...
    strain_a,s_strain_a,vabeta,creep_a);


%%  Make some plots of the result and tables

Figure6;  % slip rates
Figure7;  % block rotations
Figure8;  % residual statistics
Figure9;  % misfits per block.  If block is blank means there
Figure10; % residual velocities

% Make two tables 1) the fault slip rates and 2) the block rotation rates
% these go to the same text file.
file = 'BlockModelTables.txt';
MakeBMTable(file,slip,s_slip,omega,s_omega,creep,blocks,faults,nodes,lon,lat);


%% save the results for later

% save all the variables in a Matlab workspace
save AfterBlocks

toc