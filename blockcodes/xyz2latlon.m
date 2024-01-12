function [lat,lon,h,slat,slon,sh,covL]=xyz2latlon(x,y,z,sx,sy,sz,sxy,syz,sxz,fast)
%function
%[lat,lon,h,slat,slon,sh,covL]=xyz2latlon(x,y,z,sx,sy,sz,sxy,syz,sxz,fast)
%
% Uses the WGS84 datum (ellipsoid) to compute latitude and longitude from x,y,z data
%
% x,y,z,sx,sy,sz are in meters
% syz, syz, sxz are covariances (not correlations)
%
% fast==1 omits uncertainies and covariance for computational speed
%
% covL is covariance matrix of (lat,lon,h)
% Note that covL terms have not yet been coverted to degrees (still in
% radians).  Also covL is provided only if x,y,z are scalar.
%
% This will provide exact WGS84 coordinates if the xyz are in a WGS84
% reference frame.
%
% See WGS84 Implementation Manual from
% Eurocontrol, European Organization for the Safety of Air Navigation,
%    Brussels, Belguim and
% IfEN, Institute of Geodesy and Navigation, University FAF, Munich, Germany
% page 82.

if norm([x y z])<1000
    error('Error:  for near surface coordinates only');
end

% Earth paramters
% From Table B-1 of the above referenced document.
% Flattening of the ellipsoid
finv = 298.257223563;
f = 1./finv;
% Eccentricity of ellipsoid
e = sqrt(f*(2-f));
% semi-major axis
a =  6378137;

% requires iterative solution

% iteration 0; spherical
lat = atan((z./sqrt(x.^2 + y.^2)));
lon = atan2(y,x);

% iteration 1
nu =  a./sqrt(1-((e^2).*(sin(lat).^2)));
h = 0;
lat = atan((z./sqrt(x.^2 + y.^2)).*((1 - (e.^2).*(nu./(nu+h))).^-1));

% iteration 2 -n
tolm = .0001;  % this is the tolerance in METERS
tol = (tolm*1e-3)*pi/180/111;
lat0=lat;
dlatmax = tol*10;
itcount = 2;
while any(dlatmax > tol)
    itcount = itcount + 1;
    nu =  a./sqrt(1-e^2.*sin(lat).^2);
    h = (sqrt(x.^2 + y.^2)./cos(lat)) - nu;
    lat = atan((z./sqrt(x.^2 + y.^2)).*((1 - (e.^2)*(nu./(nu+h))).^-1));
%   dlatmax = max(abs(lat - lat0));
    dlatmax = abs(lat - lat0);
    lat0 = lat;
end

slat=[];
slon=[];
sh=[];

if ~fast
    
    slat = nan(size(lat));
    slon = nan(size(lon));
    sh = nan(size(h));
    covL = nan(length(sx),3,3);
    
    F = (1- (e.^2).*nu./(nu + h));
    
    R = sqrt(x.^2 + y.^2 + z.^2);
    A = sqrt(x.^2 + y.^2);
    
    dlon_dx = -y./(2*A.^2);
    dlon_dy =  x./(A.^2);
    
    dlat_dx = -F.*x.*z./(R.^2.*A);
    dlat_dy = -F.*y.*z./(R.^2.*A);
    dlat_dz = F.*A./(R.^2);
    
    dh_dx = x./R;
    dh_dy = y./R;
    dh_dz = z./R;
    
    for i=1:length(sx)
        
        covR = diag([sx(i).^2 sy(i).^2 sz(i).^2]);
        covR(1,2) = sxy(i);
        covR(2,1) = sxy(i);
        covR(2,3) = syz(i);
        covR(3,2) = syz(i);
        covR(1,3) = sxz(i);
        covR(3,1) = sxz(i);
        
        J = [dlat_dx(i) dlat_dy(i) dlat_dz(i);
            dlon_dx(i) dlon_dy(i) 0;
            dh_dx(i)   dh_dy(i)   dh_dz(i)];
        
        covLtemp = J*covR*J';
        
        for k=1:3
            for l=1:3
                covL(i,k,l)=covLtemp(k,l);
            end
        end
        
        slat(i,1) = sqrt(covLtemp(1,1));
        slon(i,1) = sqrt(covLtemp(2,2));
        sh(i,1) = sqrt(covLtemp(3,3));
    end
    
else
    covL=[];
end

slon = slon.*180/pi;
slat = slat.*180/pi;

lat = lat*180/pi;
lon = lon*180/pi;
