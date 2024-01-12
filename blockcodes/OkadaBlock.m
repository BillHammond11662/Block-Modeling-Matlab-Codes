function [ulon,ulat,uh]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,U,Dtop,dip,W,loneval,lateval,mode)
%[ulon,ulat,uh]=OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,U,Dtop,dip,W,loneval,lateval,mode)
%
%  nu, mu are poisson ratio and shear modulus
%  lat1,lon1,lat2,lon2 = end points of fault in degrees.
%
%  U is slip [strikeslip updip tensile] state in meters to get ulon,ulat, uh in meters
%     positive updip is thrust, negative is normal slip
%     positive strikeslip is left lateral,negative is right lateral
%
%  positive dip is down to the left while walking from p1 to p2.
%  negative dip is down to the right while walking from p1 to p2.
%
%  Dtop is distance from surface to top of fault in km
%
%  W is width of fault, i.e. distance from top to bottom along fault plane
%  in km
%
%  if mode ==0 then lateval/loneval are taken as is
%  if mode ==1 then it is assumed that desired point is near fault on left
%  side for slip rate estimation.

little=.01;
if dip<0
    lontmp = lon1;
    lattmp = lat1;
    lon1 = lon2;
    lat1 = lat2;
    lat2 = lattmp;
    lon2 = lontmp;
    dip = -dip;
    little = -little;
end

lambda = 2*mu*nu/(1-2*nu);

[delta,azstrike]=distance(lat1,lon1,lat2,lon2);
gamma = 270-azstrike;  % works for the Sumatra EQ so that fault dips down to the East.
% gamma=real(gamma);
% might be different for general cas.

L = deg2km(delta);                % in km
% WU = Mo/(1000*L*mu);            % WU in m^2
% Umag = norm(U);
% W = WU/Umag/1000;               % W in km
D = W*sin(dip*pi/180) + Dtop;   % in km
%disp(['W=' num2str(W) '; D=' num2str(D) ';L=' num2str(L)]);

dlon = loneval - lon1;
dlat = lateval - lat1;

dx = dlon*111.7.*cosd(lateval);
dy = dlat*111.7;

% this shifts the fault so that lat1,lon1 is at the origin
xoff1 = L*cosd(gamma);
yoff1 = L*sind(gamma);

% this shifts the fault so the top of the fault breaks the surface at the
% lat lons provided.
xoff2 = W*cosd(dip)*cosd(gamma+90);
yoff2 = W*cosd(dip)*sind(gamma+90);

xp = dx + xoff1 + xoff2;
yp = dy + yoff1 + yoff2;

if mode==0
    % rotate back by -gamma...
    x =  xp*cosd(gamma) + yp*sind(gamma);
    y = -xp*sind(gamma) + yp*cosd(gamma);

elseif mode==1

    x = L/2;
    y = W*cosd(dip) - little;

end

[ux,uy,uz]=DislocOkada(L,W,D,dip,U,lambda,mu,x,y,0);
uh=uz;

for i = 1:length(ux)
    ulon(i) = ux(i)*cosd(gamma) - uy(i)*sind(gamma);
    ulat(i) = ux(i)*sind(gamma) + uy(i)*cosd(gamma);
end