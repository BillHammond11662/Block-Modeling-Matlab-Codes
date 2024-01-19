function iflt=seg2faultnum(lon1,lat1,lon2,lat2,faults,nodes)
%  iflt=seg2faultnum(lon1,lat1,lon2,lat2,faults,nodes)
%

iflt=[];

P=size(faults,1);
i=1;
still=1;
while still 
    nn = faults(i,1:2);
    lonf1=nodes(nn(1),1);
    latf1=nodes(nn(1),2);
    lonf2=nodes(nn(2),1);
    latf2=nodes(nn(2),2);
    
    if ((lonf1==lon1 && latf1==lat1 && lonf2==lon2 && latf2==lat2) || ...
        (lonf1==lon2 && latf1==lat2 && lonf2==lon1 && latf2==lat1))
        iflt=i;
        still=0;
    end
    
    i=i+1;
    if i>P
        still=0;
    end
    
end