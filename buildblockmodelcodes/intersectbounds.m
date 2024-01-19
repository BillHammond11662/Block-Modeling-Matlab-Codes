function [latint,lonint]=intersectbounds(lat1,lon1,lat2,lon2,bnds)
% [lati,loni]=intersectbounds(lat1,lon1,lat2,lon2,bnds)
%
% given line segment defined by lat1,lon1,lat2,lon2
% return coordinates of point where segment intersects boundary given by 'bnds'
% returns empty values if segment does not intersect boundary
% can potentially return more than one value if segment intersects boundary
% at more than one place

%keyboard

lati=nan(4,1);
loni=nan(4,1);
ido=nan(4,1);
[ido(1),lati(1),loni(1)]=dotheyintersect(lat1,lon1,lat2,lon2,bnds(3),bnds(1),bnds(3),bnds(2));
[ido(2),lati(2),loni(2)]=dotheyintersect(lat1,lon1,lat2,lon2,bnds(3),bnds(2),bnds(4),bnds(2));
[ido(3),lati(3),loni(3)]=dotheyintersect(lat1,lon1,lat2,lon2,bnds(4),bnds(1),bnds(4),bnds(2));
[ido(4),lati(4),loni(4)]=dotheyintersect(lat1,lon1,lat2,lon2,bnds(3),bnds(1),bnds(4),bnds(1));
p=find(ido); 

if length(find(p))>1
    disp('Impossible intersection of bounds. Fix "dotheyintersect.m"');
    keyboard;
end

if isempty(p)
    latint=[];
    lonint=[];
else
    latint=lati(p);
    lonint=loni(p);
end