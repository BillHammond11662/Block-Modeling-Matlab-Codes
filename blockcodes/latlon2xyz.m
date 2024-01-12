function [x,y,z,covR]=latlon2xyz(lat,lon,h,covL)
% [x,y,z,covR]=latlon2xyz(lat,lon,h,covL);
%
%  implementation of WGS84 converstion from lat and lon to xyz.
%   lat and lon should be in degrees.
%   h is height in meters.
%
%   x,y,z are in meters in Earth fixed, Earth centered coords.
% 
%  covR is covariance matrix of x,y,z
%  covL is covariance matrix of lat,lon,h
%
%  See the WGS84 manual page 81 for the formulas
%

phi = lat*pi/180;
lambda = lon*pi/180;

finv = 298.257223563;
f = 1/finv;

% semi-major axis of the Earth
a = 6378137;

% eccentricity of the Earth
e = sqrt(f.*(2 - f));

nu = a./sqrt((1 - (e.^2).*(sin(phi)).^2));

x = (nu + h).*cos(phi).*cos(lambda);
y = (nu + h).*cos(phi).*sin(lambda);
z = (nu.*(1 - e.^2) + h).*sin(phi);

% find uncertainties

R = sqrt(x.^2 + y.^2 + z.^2);
A = sqrt(x.^2 + y.^2);

dlon_dx = -y./A.^2;
dlon_dy =  x./A.^2;
dlon_dz = zeros(length(x),1);

dlat_dx = -x.*z./(R.^2.*A);
dlat_dy = -y.*z./(R.^2.*A);
dlat_dz = A./R.^2;

dh_dx = x./R;
dh_dy = y./R;
dh_dz = z./R;

if ~isempty(covL)

    for i=1:length(x)

        J = [dlat_dx(i) dlat_dy(i) dlat_dz(i);
            dlon_dx(i) dlon_dy(i) dlon_dz(i);
            dh_dx(i)   dh_dy(i)   dh_dz(i)];

        covRtemp = J'*covL*J;

        for k=1:3
            for ll=1:3
                covR(i,k,ll)=covRtemp(k,ll);
            end;
        end;

    end;

else
    covR=[];
end;




