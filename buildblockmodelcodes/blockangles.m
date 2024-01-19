function [ angs ] = blockangles(lon,lat)
% angs  = blockangles(lon,lat)
%
% angs is a list of internal angles of the block boundary 
% assumes cartesian geometry adjusted for latitude

if (lon(1)==lon(end) && lat(1)==lat(end))
    lon(end)=[];
    lat(end)=[];
end

n=length(lon);
angs=nan(size(lon));
for i=1:n
    
    if (i>1 && i<n)
        ibef=i-1;
        imid=i;
        iaft=i+1;
    elseif i==1
        ibef=n;
        imid=1;
        iaft=2;
    elseif i==n
        ibef=n-1;
        imid=n;
        iaft=1;
    end
    
    v1 = [(lon(imid)-lon(ibef))*cosd(mean(lat));lat(imid)-lat(ibef)];
    v1 = -1*v1;
    v2 = [(lon(iaft)-lon(imid))*cosd(mean(lat));lat(iaft)-lat(imid)];
    
    angs(i) = atan2d(v1(2),v1(1)) - atan2d(v2(2),v2(1));
    
    angs(i) = mod(angs(i),360);
end


