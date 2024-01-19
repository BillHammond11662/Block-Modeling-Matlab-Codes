function iflt=seg2flt(lon1,lat1,lon2,lat2,geol)

iflt=[];
for i=1:length(geol)
    
    lonf=geol(i).lonseg;
    latf=geol(i).latseg;
    
    for j=2:length(lonf)
        
        if ((lonf(j-1)==lon1 && latf(j-1)==lat1 && lonf(j)==lon2 && latf(j)==lat2) || ...
            (lonf(j-1)==lon2 && latf(j-1)==lat2 && lonf(j)==lon1 && latf(j)==lat1))
            iflt=i;
        end
    end
    
end