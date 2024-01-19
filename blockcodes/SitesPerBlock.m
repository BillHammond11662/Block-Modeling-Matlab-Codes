function [N]=SitesPerBlock(blocks,nodes,lon,lat)
% [N]=SitesPerBlock(blocks,nodes,lon,lat)
%
% returns the number of sites in each block in the model


blocknames=fieldnames(blocks);
modelnodes=struct2cell(blocks);

for i=1:length(blocknames)

    j=modelnodes{i};
    k = find(inpolygon(lon,lat,nodes(j,1),nodes(j,2)));
    %disp([sprintf('%20s',blocknames{i}) '  '  num2str(length(k))]);  
    N(i)=length(k);
end