function [ido,lati,loni]=dotheyintersect(lat1,lon1,lat2,lon2,lat3,lon3,lat4,lon4)
%
%  ido = 1 if intersection of segments is withen ends of the segments 
%  ido = 0 otherwise
%
%  lati, loni are coordinates of intersection.  cartesian for now.
%
%  temporary code... needs to be sphericalized.

if ((lon2~=lon1 && lon4~=lon3) && (lat2~=lat1 && lat4~=lat3))
    m1 = (lat2-lat1)/(lon2-lon1);
    m2 = (lat4-lat3)/(lon4-lon3);
    
    b1= lat1 - m1*lon1;
    b2= lat3 - m2*lon3;
    
    loni = (-b1+b2)/(m1-m2);
    lati = m1*loni + b1;
elseif lon2==lon1
    loni=lon2;
    lati = lat3 + (lon2-lon3)*(lat4 - lat3)/(lon4 - lon3);
elseif lon3==lon4
    loni = lon3;
    lati = lat1 + (lon3-lon1)*(lat2 - lat1)/(lon2 - lon1);
elseif lat2==lat1
    lati=lat2;
    loni = lon3 + (lat3 - lat2)*(lon4 - lon3)/(lat4 - lat3);
elseif lat4==lat3
    lati=lat4;
    loni = lon1 + (lat4 - lat1)*(lon2 - lon1)/(lat2 - lat1);
end

% see if there is intersection
ido=0;
if (min([lon1 lon2]) <= loni && loni <= max([lon1 lon2]) && ...
     min([lon3 lon4]) <= loni && loni <= max([lon3 lon4]) && ... 
     min([lat1 lat2]) <= lati && lati <= max([lat1 lat2]) && ...
      min([lat3 lat4]) <= lati && lati <= max([lat3 lat4]))
    ido=1;
end

% see if intersection is at endpoint (i.e. pt1=pt3 or pt1=pt4 or pt2==pt3
% or pt2==pt4
% if ido
%     if ((lon1==lon3 && lat1==lat3) || ...
%         (lon1==lon4 && lat1==lat4) || ...
%         (lon2==lon3 && lat2==lat3) || ...
%         (lon2==lon4 && lat2==lat4))
%         ido = 0;
%     end
% end
