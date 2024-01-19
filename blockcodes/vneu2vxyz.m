function [vx,svx,vy,svy,vz,svz,covxyz]=vneu2vxyz(vn,svn,ve,sve,vu,svu,covneu,lat,lon)
% function [vx,svx,vy,svy,vz,svz,covxyz]=vneu2vxyz(vn,svn,ve,sve,vu,svu,covneu,lat,lon)
%
% make sure lat, lon are in degrees;
%
% covneu is in format (covne, coveu, covnu) in columns
%
% note: obliquity not taken into account.

% to make up velocity negative (as per p. 155 "Plate Tectonics: How it Works" Cox & Hart)
vu = -vu;

lat = lat*pi/180;
lon = lon*pi/180;

Tnx = -sin(lat).*cos(lon);
Tex = -sin(lon);
Tdx = -cos(lat).*cos(lon);
Tny = -sin(lat).*sin(lon);
Tey =  cos(lon);
Tdy = -cos(lat).*sin(lon);
Tnz =  cos(lat);
Tez =  zeros(size(lat));
Tdz = -sin(lat);

vx=nan(size(lat));
vy=nan(size(lat));
vz=nan(size(lat));
svx=nan(size(lat));
svy=nan(size(lat));
svz=nan(size(lat));

for i =1:length(lat)
    T = [Tnx(i) Tny(i) Tnz(i);
        Tex(i) Tey(i) Tez(i);
        Tdx(i) Tdy(i) Tdz(i)];

    Tinv = inv(T);

    VNEU = [vn(i) ve(i) vu(i)]';
    VXYZ = Tinv*VNEU;

    COVNEU = diag([svn(i);sve(i);svu(i)].^2);

    if ~isempty(covneu)
        COVNEU(1,2) = covneu(i,1);
        COVNEU(2,1) = covneu(i,1);
        COVNEU(2,3) = covneu(i,2);
        COVNEU(3,2) = covneu(i,2);
        COVNEU(1,3) = covneu(i,3);
        COVNEU(3,1) = covneu(i,3);
    end

    covxyz = Tinv*COVNEU*Tinv';

    vx(i) = VXYZ(1);
    vy(i) = VXYZ(2);
    vz(i) = VXYZ(3);
    %
    % 	svx(i) = sqrt((Tinv(1,1)*svn(i)).^2 + (Tinv(1,2).*sve(i)).^2 + (Tinv(1,3).*svu(i)).^2);
    % 	svy(i) = sqrt((Tinv(2,1)*svn(i)).^2 + (Tinv(2,2).*sve(i)).^2 + (Tinv(2,3).*svu(i)).^2);
    % 	svz(i) = sqrt((Tinv(3,1)*svn(i)).^2 + (Tinv(3,2).*sve(i)).^2 + (Tinv(3,3).*svu(i)).^2);

    svx(i) = sqrt(covxyz(1,1));
    svy(i) = sqrt(covxyz(2,2));
    svz(i) = sqrt(covxyz(3,3));
    
end;
