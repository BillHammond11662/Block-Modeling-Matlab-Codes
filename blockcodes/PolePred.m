function [vn,ve]=PolePred(E,lat,lon)
% function [vn,ve]=PolePred(E,lat,lon)
%
% lat and lon are in DEGREES
% E is Euler vector where
% E = [OMEGA*cos(PHI_p)*cos((pi/2)-THETA_p);
%      OMEGA*sin(PHI_p)*cos((pi/2)-THETA_p);
%      OMEGA*sin((pi/2)-THETA_p)];
% where OMEGA is rotation rate, PHI is longitude of pole, THETA is CO-latitude
% of pole

[n,m]=size(lon);
ve = NaN(n,m);
vn = NaN(n,m);

for i=1:n
    for j=1:m

        lonr = lon(i,j).*pi/180;
        latr = lat(i,j).*pi/180;
        R = 6371009.3;

        % find the velocity at the site predicted by the pole.
        r = [R.*cos(lonr).*cos(latr);
            R.*sin(lonr).*cos(latr);
            R.*sin(latr)];
        vpred = cross(E,r);

        % convert to vn,ve,vu
        T = [-sin(latr)*cos(lonr)  -sin(latr)*sin(lonr)  cos(latr);
            -sin(lonr)                 cos(lonr)                0;
            -cos(latr)*cos(lonr)  -cos(latr)*sin(lonr)  -sin(latr)];

        Vneu = T*vpred;
        vn(i,j)= Vneu(1);
        ve(i,j)= Vneu(2);

    end
end