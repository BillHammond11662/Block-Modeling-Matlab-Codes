function j=FindBlock(blocks,nodes,lon,lat)
% function j=FindBlock(blocks,nodes,lon,lat)
%
%  gives block number from a particular lat,lon
blocknames = fieldnames(blocks);
M = length(blocknames);
j=[];
for jj=1:M
    eval(['jb=blocks.' char(blocknames(jj)) ';']);
    xb = nodes(jb,1);
    yb = nodes(jb,2);
    if find(inpolygon(lon,lat,xb,yb))
        j=[j;jj];
    end
end;