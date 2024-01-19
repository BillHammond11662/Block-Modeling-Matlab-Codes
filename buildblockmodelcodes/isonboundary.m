function ison = isonboundary(lat1,lon1,lat2,lon2,bnds)
% ison = isonboundary(lat1,lon1,lat2,lon2,bnds)
%

ison=0;

if lat1==lat2 && (lat1==bnds(3) || lat1==bnds(4))
    ison=1;
elseif lon1==lon2 && (lon1==bnds(1) || lon1==bnds(2))
    ison=1;
end