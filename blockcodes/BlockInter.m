function [omega,s_omega,slip,s_slip,strain,s_strain,covm,vnpred,vepred,chi2,dof,ilist,mse]=...
    BlockInter(lon,lat,ve,vn,vu,se,sn,su,nodes,blocks,faults,nu,mu,alpha,beta,gamma,zeta,L,...
    omega_a,s_omega_a,slip_a,s_slip_a,strain_a,s_strain_a,vabeta,creep)
% [omega,s_omega,slip,s_slip,strain,s_strain,covm,vnpred,vepred,chi2,dof,ilist,mse]=...
%     BlockInter(lon,lat,ve,vn,vu,se,sn,su,nodes,blocks,faults,nu,mu,alpha,beta,gamma,zeta,L,...
%     omega_a,s_omega_a,slip_a,s_slip_a,strain_a,s_strain_a,vabeta,creep)
%
%  INPUTS
%  lon,lat == coordinates of vel data
%  ve,vn,vu,se,sn,su == vel data and uncertainty in m/yr
%  nodes is lise of lon,lat coordates of each block boundary defining node
%  blocks defines the block boundaryies from the nodes
%  faults defines the faults which will be constrained to slip
%
%  if vu is NaNs then no vertical constraints are applied
%
%  alpha = slip rate consistency constraint (m/yr should be very small, try .00001)
%  beta = a priori uncertainty on rotation vector (rad/yr)
%  gamma = a priori uncertainty on slip rates (m/yr)
%  zeta = a priori uncertainty on strain rates (strains/yr)
%
%  L = maximum number of nearnest fault segments.
%
%  omega_a,s_omega_a   = a priori omega value, a priori uncertainty in omega
%  slip_a,s_slip_a     = a priori slip rate value, a priori uncertainty in slip rate
%  strain_a,s_strain_a = a priori strain rate value, a priori uncertainty in strain rate
%
%  Each of the above must have size of omega,slip,strain output matrices
%  and have NaNs in all entries that are not provided.
%  Also note,  for omega_a and strain_a you need to provide all three
%  parameters. However, slip_a can take just one or all two.
%
%  creep is a K by 2 matrix, where K is the number of faults whose creep
%  rate is imposed.  First column is the fault number, second column is the
%  creep rate in m/yr, positive left lateral.
%
%  OUTPUTS:
%  omega is a list of rotation vectors for each block
%  s_omega is uncertainty in omega
%
%  slip is a list of slip vectors, one for each fault
%  s_slip is uncertainty in slip
%  slip(:,1) are strike slip rates, slip(:,2) are dip slip rates (in along dip direction not projected to horizontal)
%  positive slip(:,1) is sinistral, positive slip(:,2) is thrust
%  negative slip(:,1) is normal, negative slip(:,2) is normal
%
%  strain is strain in each block (NaNs if zeta = 0); e_phiphi e_thetaphi e_thetatheta
%  s_strain is uncertainty in strains
%
%  covm is covariance matrix with dimensions (3*M + 2*P) by (3*M + 2*P)
%  where M is the number of blocks and
%        P is the number of faults
%
%  positive fault dip is down to the left while walking from p1 to p2.
%  negative fault dip is down to the right while walking from p1 to p2.
%
%  to mimic the case where there is no strain accumulation, set the elastic
%  layer thickness in "faults" to a very small number.
%
%  vabeta is a flag that changes the damping constraint on block rotation
%  to damping on the vertical axis component (spin) of the rotation. 
%       vabeta = 1 is damping on vertical axis rotation
%       otherwise damping is on omega
%
%  chi2 is misfit.
%  dof is degrees of freedom (number of data minue number of model parameters)
%
%  ilist are the indexes of lon,lat,ve,vn,se,sn of GPS stations that line
%  on some block and were used as data in the inversion.
%  
%  Note vu, su are not used at this time.

r0 = 6.378e6;

blocknames = fieldnames(blocks);
M = length(blocknames);
[P,~] = size(faults);

scal = 1e7;  % scales the slip rate parameters for matrix balancing
% scal = 1;
gamma = gamma/scal;  % adjust gamma since the referred to parameters are scaled.

iln=find(lon>180);
lon(iln)=lon(iln)-360;

iln2=find(nodes(:,1)>180);
nodes(iln2,1)=nodes(iln2,1)-360;

% from fault data, determine the fault midpoints (in degrees)
plat=nan(P,1);
plon=nan(P,1);
for i=1:P
    [rng(i),az(i)]=distance(nodes(faults(i,1),2),nodes(faults(i,1),1), ...
        nodes(faults(i,2),2),nodes(faults(i,2),1));
    [plat(i,1),plon(i,1)]=reckon(nodes(faults(i,1),2),nodes(faults(i,1),1),rng(i)/2,az(i));
end

% xyz coordinates of fault midpoints for later use
px=NaN(P,1);
py=NaN(P,1);
pz=NaN(P,1);
for k=1:P
    [px(k,1),py(k,1),pz(k,1)]=latlon2xyz(plat(k),plon(k),0,zeros(3,3));
end
p = [px py pz];

% make a list of block centroids
b=zeros(M,2);
blon=nan(M,1);
blat=nan(M,1);
for i=1:M
    jb=blocks.(blocknames{i});
    xb = nodes(jb,1);
    yb = nodes(jb,2);
    [blon(i,1),blat(i,1)]=centroid(polyshape(xb,yb,'KeepCollinearPoints',true));
    [b(i,1),b(i,2),b(i,3)]=latlon2xyz(mean(yb(1:(end-1))),mean(xb(1:(end-1))),0,zeros(3,3));
end

%% Define the velocity data constraint

disp('Adding velocity data constraint on block rotation...');

rz=[0;0;1];
d = [];
G = [];
ilist = [];
vsmall = 1e-6;

for i=1:length(lon)
    
    j=FindBlock(blocks,nodes,lon(i),lat(i));
    
    while length(j)>1
        %disp(['Warning: Site ' num2str(i) ' is on more than one block. Choosing one.']);
        lon(i) = lon(i)+vsmall;
        lat(i) = lat(i)+vsmall;
        j=FindBlock(blocks,nodes,lon(i),lat(i));
    end
    
    if ~isempty(j)  % the velocity has to be inside some block
        
        [x,y,z]=latlon2xyz(lat(i),lon(i),0,zeros(3,3));
        r=[x;y;z];
        
        % unit vectors in north & east direction
        en = cross(cross(r,rz),r);
        en = en/norm(en);
        ee = cross(en,r);
        ee = ee/norm(ee);
        gn=cross(-en,r)';
        ge=cross(-ee,r)';
        
        ilist = [ilist;i];   % the list of sites that found a block
        d = [d;
            vn(i)
            ve(i)];
        
        G = [G;
            zeros(1,3*(j-1))   gn  zeros(1,3*(M-j));
            zeros(1,3*(j-1))   ge  zeros(1,3*(M-j))];
    end
    
end

N = length(ilist);

% G is now 2N x 3M

%% Add the block secular strain effects

if ~(zeta==0 || isnan(zeta) || isempty(zeta))
    
    disp('Adding block secular strain parameters...');
    
    %G2 = [];  
    G2=zeros(2*length(ilist),3*M);
    for i=1:length(ilist)
        
        [x,y,z]=latlon2xyz(lat(ilist(i)),lon(ilist(i)),0,zeros(3,3));
        r=[x;y;z];
        en = cross(cross(r,rz),r);
        en = en/norm(en);% unit vector in north direction
        ee = cross(en,r);
        ee = ee/norm(ee);% unit vector in east direction
        
        j=FindBlock(blocks,nodes,lon(ilist(i)),lat(ilist(i)));
        
        theta0 = 90-lat(ilist(i));
        delphi = (lon(ilist(i)) - blon(j))*pi/180;
        deltheta = -(lat(ilist(i)) - blat(j))*pi/180;
        
        G2(2*(i-1)+1,3*(j-1)+1:3*(j-1)+3)= [ 0                      -r0*sind(theta0)*delphi -r0*deltheta ];
        G2(2*(i-1)+2,3*(j-1)+1:3*(j-1)+3)= [ r0*sind(theta0)*delphi  r0*deltheta             0           ];
    end
    G = [G G2];

end

% if zeta non zero or blank
% G is now 2N x (3M + 3M)
% else
% G is still 2N x 3M


%% Add the fault-related elastic strain accumulation effect

G = [G zeros(2*N,2*P)];

Uss = [1 0 0];
Un  = [0 1 0];  % note this is actually thrust but that's ok
Dtop = 0;

disp('Adding elastic strain accumulation effect...');

L = min([L size(faults,1)]);

if ~isempty(faults)
    for i=1:N
        dist=sqrt((lon(ilist(i))-plon).^2 + (lat(ilist(i))-plat).^2);
        [~,isrt]=sort(dist);
        
        for kk=1:L
            knum = isrt(kk);
            
            if zeta~=0
                k1 = 6*M + 2*(knum-1) + 1;
            else
                k1 = 3*M + 2*(knum-1) + 1;
            end
            
            lat1 = nodes(faults(knum,1),2);
            lon1 = nodes(faults(knum,1),1);
            lat2 = nodes(faults(knum,2),2);
            lon2 = nodes(faults(knum,2),1);
            dip  = faults(knum,4);
            W    = abs(faults(knum,3)/sind(dip));
            
            [glonss,glatss,~]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Uss,Dtop,dip,W,lon(ilist(i)),lat(ilist(i)),0);
            [glonn,glatn,~]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Un,Dtop,dip,W,lon(ilist(i)),lat(ilist(i)),0);
            
            G(2*(i-1)+1,    k1) = -glatss*scal;
            G(2*(i-1)+1,  k1+1) = -glatn*scal;
            G(2*(i-1)+2,    k1) = -glonss*scal;
            G(2*(i-1)+2,  k1+1) = -glonn*scal;
        end
    end
end

% if zeta non zero or blank
% G is now 2N x (3M + 3M + 2P)
% else
% G is still 2N x (3M + 2P)

%% Add the slip rate constraint

disp('Adding slip rate/block rotation consistency constraint...');

if zeta~=0
    G = [G;zeros(2*P,6*M+2*P)];
else
    G = [G;zeros(2*P,3*M+2*P)];
end
d = [d;zeros(2*P,1)];

% dcord = .5;  % in km

for k=1:P
    
    % find the blocks involved
    
    [j]=Fault2Blocks(nodes,blocks,faults,k);
    
    %    evaluate just off the fault.
    lat1 = nodes(faults(k,1),2);
    lon1 = nodes(faults(k,1),1);
    lat2 = nodes(faults(k,2),2);
    lon2 = nodes(faults(k,2),1);
    dip  = faults(k,4);
    W    = abs(faults(k,3)/sind(dip));
    
    [rngk,~]=distance(lat1,lon1,lat2,lon2);
    dcord = min([deg2km(rngk)/100 .5]);
    
    % find the positions slightly off the fault on both blocks
    % to get relative slip rate, i.e motion of block 1 wrt block 2
    % determine which block is block 2 and which is 1
    [lateval, loneval]=reckon(plat(k),plon(k),km2deg(dcord),az(k)-90);
    jtest=FindBlock(blocks,nodes,loneval,lateval);
    if isempty(jtest)
        [lateval, loneval]=reckon(plat(k),plon(k),km2deg(dcord),az(k)+90);
        jtest=FindBlock(blocks,nodes,loneval,lateval);
    end
    if jtest==j(1)
        loneval1=loneval;
        lateval1=lateval;
        [lateval2, loneval2]=reckon(plat(k),plon(k),km2deg(dcord),az(k)+90);
    elseif jtest==j(2)
        loneval2=loneval;
        lateval2=lateval;
        [lateval1, loneval1]=reckon(plat(k),plon(k),km2deg(dcord),az(k)+90);
    else
        disp('Could not figure out which is block 1 or 2. Aborting this model.');
        omega=[]; s_omega=[];
        slip=[]; s_slip=[];
        strain=[]; s_strain=[];
        covm=[];
        vnpred=[]; vepred=[];
        chi2=[];dof=[];
        ilist=[];mse=[];
        return;
        %         omega=[];s_omega,slip,s_slip,strain,s_strain,covm,vnpred,vepred,chi2,dof,ilist,mse
        %         return;
    end
        
    % try to do this with option 1
    Gss=[];Gn=[];

    [glonss,glatss,ghss]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Uss,Dtop,dip,W,loneval1,lateval1,0);
    [glonn,glatn,ghn]   =OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Un,Dtop,dip,W,loneval1,lateval1,0);
    
    glonss=real(glonss);
    glatss=real(glatss);
    ghss=real(ghss);
    ghn=real(ghn);
    
    [Gss(1,1),~,Gss(2,1),~,Gss(3,1),~]=...
        vneu2vxyz(glatss,1,glonss,1,ghss,1,[],plat(k),plon(k));
    [Gn(1,1),~,Gn(2,1),~,Gn(3,1),~]=...
        vneu2vxyz(glatn,1,glonn,1,ghn,1,[],plat(k),plon(k));
    Gss1=real(Gss);
    Gn1=real(Gn);
    
    Gss=[];Gn=[];
    [glonss,glatss,ghss]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Uss,Dtop,dip,W,loneval2,lateval2,0);
    [glonn,glatn,ghn]   =OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Un,Dtop,dip,W,loneval2,lateval2,0);
    
    glonss=real(glonss);
    glatss=real(glatss);
    ghss=real(ghss);
    ghn=real(ghn);
    
    [Gss(1,1),~,Gss(2,1),~,Gss(3,1),~]=...
        vneu2vxyz(glatss,1,glonss,1,ghss,1,[],plat(k),plon(k));
    [Gn(1,1),~,Gn(2,1),~,Gn(3,1),~]=...
        vneu2vxyz(glatn,1,glonn,1,ghn,1,[],plat(k),plon(k));
    Gss2=real(Gss);
    Gn2=real(Gn);
    
    dGss = Gss1 - Gss2; 
    dGn = Gn1 - Gn2;     
        
    if (norm(dGss)<.8)
        disp(['Norm of dGss<0.8 : Fault #' num2str(k)]);
    end
    if (norm(dGn)<.8)
        disp(['Norm of dGn<0.8 : Fault #' num2str(k)]);
    end
    
    pz=[0;0;1];
    en = cross(cross(p(k,:)',pz),p(k,:)');
    en = en/norm(en);% unit vector in north direction
    ee = cross(en,p(k,:)');
    ee = ee/norm(ee);% unit vector in east direction
       
    gn1=cross(-en,p(k,:)')';
    gn2=cross(en,p(k,:)')';
    ge1=cross(-ee,p(k,:)')';
    ge2=cross(ee,p(k,:)')';
    
    ij1 = 3*(j(1)-1)+1:3*(j(1)-1)+3;
    ij2 = 3*(j(2)-1)+1:3*(j(2)-1)+3;
    ik = 2*(k-1)+1;
    
    G(2*N+ik,  ij1)= gn1';
    G(2*N+ik,  ij2)= gn2';
    G(2*N+ik+1,ij1)= ge1';
    G(2*N+ik+1,ij2)= ge2';
    
    theta0 = 90-plat(k);
    delphi1 = (plon(k) - blon(j(1)))*pi/180;
    delphi2 = (plon(k) - blon(j(2)))*pi/180;
    deltheta1 = (-plat(k) + blat(j(1)))*pi/180;
    deltheta2 = (-plat(k) + blat(j(2)))*pi/180;

    jc=[];
    if ~isempty(creep)
       jc=find(k==creep(:,1));
    end
 
    if zeta~=0
        
        G(2*N+ik,  3*M+ij1)= [0                        -r0*sind(theta0)*delphi1 -r0*deltheta1];
        %G(2*N+ik,  3*M+ij2)= [0                         r0*sind(theta0)*delphi2  r0*deltheta2];
        G(2*N+ik,  3*M+ij2)= [0                         -r0*sind(theta0)*delphi2  -r0*deltheta2];

        G(2*N+ik+1,3*M+ij1)= [-r0*sind(theta0)*delphi1 -r0*deltheta1             0];
%         G(2*N+ik+1,3*M+ij2)= [ r0*sind(theta0)*delphi2  r0*deltheta2             0];
        G(2*N+ik+1,3*M+ij2)= [ -r0*sind(theta0)*delphi2  -r0*deltheta2             0];

        
        G(2*N+ik,  6*M+ik)  = -dot(dGss,en)*scal;
        G(2*N+ik,  6*M+ik+1)= -dot(dGn,en)*scal;
        G(2*N+ik+1,6*M+ik)  = -dot(dGss,ee)*scal;
        G(2*N+ik+1,6*M+ik+1)= -dot(dGn,ee)*scal;
        
        if ~isempty(jc)
            
            [vxkap,~,vykap,~,vzkap,~,~] = vneu2vxyz(cosd(az(k)),1,sind(az(k)),1,0,1,[],plat(k),plon(k));
            vkappa = mean(creep(jc,2))*[vxkap vykap vzkap]'; % this is the xyz creep vector
            
            d(2*N+ik)= -dot(vkappa,en);  % should be north component of xyz vector kappa            
            d(2*N+ik+1)= -dot(vkappa,ee);  % should be east component of xyz vector
            
        end
    else
        
        G(2*N+ik,  3*M+ik)  = -dot(dGss,en)*scal;
        G(2*N+ik,  3*M+ik+1)= -dot(dGn,en)*scal;
        G(2*N+ik+1,3*M+ik)  = -dot(dGss,ee)*scal;
        G(2*N+ik+1,3*M+ik+1)= -dot(dGn,ee)*scal;
        
        if ~isempty(jc)
            
            [vxkap,~,vykap,~,vzkap,~,~] = vneu2vxyz(cosd(az(k)),1,sind(az(k)),1,0,1,[],plat(k),plon(k));
            vkappa = mean(creep(jc,2))*[vxkap vykap vzkap]'; % this is the xyz creep vector
            
            d(2*N+ik)= -dot(vkappa,en);  % should be north component of xyz vector kappa            
            d(2*N+ik+1)= -dot(vkappa,ee);  % should be east component of xyz vector
            
        end
    end
    
end

% if zeta non zero or blank
% G is now (2N + 2P) x (6M + 2P)
% else
% G is still (2N +2P) x (3M + 2P)


%%   Add apriori uncertainty - set the weighting 

disp('Adding weights on existing equations using alpha, beta, gamma, zeta... (data vector focussed regularization)');

clear W;
W=reshape([1./se(ilist).^2 1./sn(ilist).^2]',1,2*length(ilist))';

if (~isnan(zeta) && zeta~=0)
    W = [W;...                          % data weighting
        (1./alpha.^2)*ones(2*P,1);...   % slip rate consistency weighting
        (1./beta.^2)*ones(3*M,1);...    % rotation rate weighting
        (1./zeta.^2)*ones(3*M,1);...    % strain rate weighting
        (1./gamma.^2)*ones(2*P,1)];     % slip rate weighting
else
    W = [W;...                          % data weighting
        (1./alpha.^2)*ones(2*P,1);...   % slip rate consistency weighting
        (1./beta.^2)*ones(3*M,1);...    % rotation rate weighting
        (1./gamma.^2)*ones(2*P,1)];     % slip rate weighting.  
end


%% apply apriori constraints on model parameters

disp('Adding damping constraints ... (model parameter focussed regularization)');
%disp('These push model parameters towards certain values (e.g., zero or geologic rate');

[N2,M2]=size(G);
Gadd=eye(M2);

G=[G;Gadd];
d=[d;zeros(M2,1)];

for i=1:M
    if ~isnan(omega_a(i,1))
        j = N2+3*(i-1);
        d(j+1)=omega_a(i,1);
        d(j+2)=omega_a(i,2);
        d(j+3)=omega_a(i,3);

        W(j+1)=(1./s_omega_a(i,1)).^2;
        W(j+2)=(1./s_omega_a(i,2)).^2;
        W(j+3)=(1./s_omega_a(i,3)).^2;
    end
end

if (~isnan(zeta) && zeta~=0)
    for i=1:M
        if ~isnan(strain_a(i,1))
            j = N2+3*M + 3*(i-1);
            d(j+1)=strain_a(i,1);
            d(j+2)=strain_a(i,2);
            d(j+3)=strain_a(i,3);

            W(j+1)=(1./s_strain_a(i,1)).^2;
            W(j+2)=(1./s_strain_a(i,2)).^2;
            W(j+3)=(1./s_strain_a(i,3)).^2;
        end
    end
end

for i=1:P
    
    if (~isnan(zeta) && zeta~=0)
       j = N2+6*M+2*(i-1);
    else
       j = N2+3*M+2*(i-1);
    end
    if ~isnan(slip_a(i,1)) && ~isnan(s_slip_a(i,1))
        d(j+1)=slip_a(i,1)/scal;
        W(j+1)=1./(s_slip_a(i,1)/scal).^2;  % bill 2022-05-04
    end
    if ~isnan(slip_a(i,2)) && ~isnan(s_slip_a(i,2))
        d(j+2)=slip_a(i,2)/scal;
        W(j+2)=1./(s_slip_a(i,2)/scal).^2; % bill 2022-05-04
    end
end

% if zeta non zero or blank
% G is now (2N + 2P + 6M + 2P) x (6M + 2P)
% else
% G is still (2N + 2P + 3M + 2P) x (3M + 2P)

if vabeta==1
    
    % 1- take out 3M rows owing for damping rotation, EXCEPT for rows that
    % have omega_a not nan
    
    % 2- add M rows for wx*rx + wy*ry + wz*rz = 0, EXCEPT for rows
    %     that have omega_a not nan
    
    %1
    bout=find(isnan(omega_a(:,1)) & isnan(omega_a(:,2)) & isnan(omega_a(:,3)));
    iout = reshape([3*(bout-1)+1 3*(bout-1)+2 3*(bout-1)+3]',3*length(bout),1);
    G(N2+iout,:)=[];
    W(N2+iout)=[];
    d(N2+iout)=[];
    
    %J=length(iout);
    
% if zeta non zero or blank
% G is now (2N + 2P + 6M + 2P - J) x (6M + 2P)
% else
% G is still (2N + 2P + 3M + 2P - J) x (3M + 2P)   

    %2
    for i=1:M
        
        if ~isempty(intersect(i,bout))
            %ndlist=getfield(blocks,char(blocknames(i)));
%             ndlist=blocks.(char(blocknames(i)));
%             lon0=mean(nodes(ndlist,1));
%             lat0=mean(nodes(ndlist,2));
%             
            [rx,ry,rz,~]=latlon2xyz(blat(i),blon(i),0,[]);
            rs = [rx ry rz]./norm([rx ry rz]);
            rbeg = zeros(1,3*(i-1));
            rend = zeros(1,3*(M-i));
            
            G=[G;rbeg rs rend zeros(1,M2-3*M)];
            W=[W; 1./beta.^2 ];
            d=[d; 0];
        end
        
    end

% if zeta non zero or blank
% G is now (2N + 2P + 6M + 2P - J + J/3) x (6M + 2P)
% else
% G is still (2N + 2P + 3M + 2P - J + J/3) x (3M + 2P)  
    
end

%% Compute the inverse

% keyboard;

disp('Solving...');

[nn,mm]=size(G);
[ww,vv]=size(W);
[dd,~]=size(d);
disp(['G matrix is ' num2str(nn) ' x ' num2str(mm)]);
disp(['W matrix is ' num2str(ww) ' x ' num2str(vv)]);
disp(['d vector is ' num2str(dd) ' x 1']);

if ww~=nn
    disp('Problem with matrix equation weights. ww~=nn');
    omega=[]; s_omega=[];
    slip=[]; s_slip=[];
    strain=[]; s_strain=[];
    covm=[];
    vnpred=[]; vepred=[];
    chi2=[];dof=[];
    ilist=[];mse=[];
    return;
end
if any(isnan(W))
    disp('Problem: NaNs in weight vector. Aborting this inversion.');
    omega=[]; s_omega=[];
    slip=[]; s_slip=[];
    strain=[]; s_strain=[];
    covm=[];
    vnpred=[]; vepred=[];
    chi2=[];dof=[];
    ilist=[];mse=[];
    return;
end

% W=diag(W);
% GWG = G'*W*G;
% Minv = inv(GWG)*G'*W;
% m = Minv*d;
% covm = Minv*covd*Minv';
% sm = sqrt(diag(covm));

% use this
% [m,sm,~,covm] = lscov(G,d,W);
% mse=[];

% or this 
mse=[];
Gg = (G'*diag(W)*G)\(G'*diag(W));
m = Gg*d;
covd = diag(1./W);
covm = Gg*covd*Gg';
sm = sqrt(diag(covm));

disp(['Gg matrix is ' num2str(size(Gg,1)) ' x ' num2str(size(Gg,2))]);
disp(['m vector is ' num2str(size(m,1)) ' x ' num2str(size(m,2))]);
disp(['sm vector is ' num2str(size(sm,1)) ' x ' num2str(size(sm,2))]);


%% Find the model predictions and misfit

dpred=G*m;
dpred2=reshape(dpred(1:2*N)',2,N)';
vepred = nan(size(ve));
vnpred = nan(size(vn));
vnpred(ilist) = dpred2(:,1);
vepred(ilist) = dpred2(:,2);

if zeta~=0
    m(6*M+1:end)=m(6*M+1:end)*scal;
else
    m(3*M+1:end)=m(3*M+1:end)*scal;
end

% covd = diag(1./diag(W));

for i=1:M
    jj = 3*(i-1)+1;
    omega(i,1:3) = m(jj:jj+2)';
    s_omega(i,1:3) = sm(jj:jj+2)';

    if zeta~=0
        strain(i,1:3)=m(3*M+jj:3*M+jj+2)';
        s_strain(i,1:3)=sm(3*M+jj:3*M+jj+2)';
    else
        strain(i,1:3)=NaN;
        s_strain(i,1:3)=NaN;
    end
end

for k=1:P
    if zeta~=0
        ik = 6*M+2*(k-1)+1;
    else
        ik = 3*M+2*(k-1)+1;
    end
    slip(k,1) = m(ik);
    s_slip(k,1) = sm(ik)*scal;
    slip(k,2) = m(ik+1);
    s_slip(k,2) = sm(ik+1)*scal;
end

% dof = 2*N - 6*M;
%[~,M3]=size(G);
%dof = N3 - M3;
dof = 2*length(ilist) - 3*M;
sd = reshape([sn(ilist) se(ilist)]',2*N,1);
chi2 = sum(((d(1:2*N) - dpred(1:2*N))./sd).^2);

% check out ok

%chi2 = sum(  [((ve(ilist) - vepred(ilist))./se(ilist)).^2 ;  ((vn(ilist) - vnpred(ilist))./sn(ilist)).^2]    )



