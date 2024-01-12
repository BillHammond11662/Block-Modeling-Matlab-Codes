% Run a Block Model

clear;
close all;
tic

addpath('../blockcodes');
addpath('../blockcodes/figurecodes');


%%  Get model geometries. Blocks and faults. Make some plots.

disp('Reading in block model from the .mat file.')
load blockmodel.mat;
% the matlab workspace file blockmodel.mat has variables 'nodes', 'blocks',
% and 'faults'.  These are the minimum necessary to define the block model
% geometries. 

% Uncomment/comment out the next line to run or not run the code SimplifyBlocks.m 
% that reduces the number of blocks from 62 to 24.
% SimplifyBlocks;

disp('Performing Model Checks...')
CheckNodes(nodes, faults, blocks);
% Note that "CheckNodes" gives a warning that node 98 is very close to node 99.  
% This is OK in this case so we'll forge ahead

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


%% Run the block model, i.e., solve for block rotation and slip rates.

blocknames = fieldnames(blocks);
M=size(blocknames,1);
P=size(faults,1);

% specify some basic parameters for the block model
nu = .25;           % Poisson's ration of crust (unitless)      
mu = 3e10;          % shear modulus of crust (Pa)

% these regularize the model by specifying weights or damping for each model parameter
alpha = 0.00001;    % apriori uncertainty in slip rate consistency constraint, units in "m/yr"
beta = 3e-9;        % apriori uncertainty in omega in units "rad/yr" for weighting
gamma = .0020;       % apriori uncertainty in slip rates in units "m/yr" for weighting
vabeta = 1;         % if vabeta==1 then vertical axis spin is damped instead of total rotation Euler pole
L=6;                % L is max number of faults providing interseismic strain from locking.

% In this case we will not account for active creep on faults
creep=[];
creep_a=[];

% In this case apriori values of slip rates are not applied (so we put
% nans).  Can do this to try, e.g., a geologic rate
slip_a=nan(P,2);
s_slip_a=nan(P,2);

% slip_a(131:137,1)=-.002;
% s_slip_a(131:137,1)=.0002;
% slip_a(131:137,2)=0;
% s_slip_a(131:137,2)=.0002;

% In this case apriori values of strain rates are not applied (so we put
% nans)
zeta = 0;
strain_a = nan(M,3);
s_strain_a = nan(M,3);

% In this case apriori values of block rotation rates are not applied (so we put
% nans)
omega_a = nan(M,3);
s_omega_a = nan(M,3);

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
Figure9;  % misfits per block
Figure10; % residual velocities

% Make two tables 1) the fault slip rates and 2) the block rotation rates
% these go to the same text file.
file = 'BlockModelTables.txt';
MakeBMTable(file,slip,s_slip,omega,s_omega,creep,blocks,faults,nodes,lon,lat);


%% save the results for later

save AfterBlocks

toc