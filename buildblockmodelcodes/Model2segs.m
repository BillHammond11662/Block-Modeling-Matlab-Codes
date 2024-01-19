function [lonlatseg,ndseg]=Model2segs(blocks,nodes)
% [lonlatseg,ndseg]=Model2segs(blocks,nodes)
% 
% provides complete list of block boundary segments in block model
% each row of lonlatseg is a segment lon/lats: [lon1 lat1 lon2 lat2]
% 
% each row of ndseg are the node numbers and has same number of rows as lonlatseg
%

% create list of segments
bnames=fieldnames(blocks);
M=length(bnames);
lonlatseg=[];
ndseg=[];
for i=1:M
    
    bname = char(bnames(i));
    eval(['nn=blocks.' bname ';']);
    
    for j=1:(length(nn)-1)
        
        lon1=nodes(nn(j),1);
        lat1=nodes(nn(j),2);
        lon2=nodes(nn(j+1),1);
        lat2=nodes(nn(j+1),2);
        
        %only add the segment if it is not already present in the list
        if i==1
            lonlatseg=[lonlatseg;lon1 lat1 lon2 lat2]; %#ok<*AGROW>
            ndseg=[ndseg;nn(j) nn(j+1)]; %#ok<*AGROW>
        else
            isit = find((lonlatseg(:,1)==lon1 & lonlatseg(:,2)==lat1 & lonlatseg(:,3)==lon2 & lonlatseg(:,4)==lat2 ) | ...
                        (lonlatseg(:,3)==lon1 & lonlatseg(:,4)==lat1 & lonlatseg(:,1)==lon2 & lonlatseg(:,2)==lat2 ), 1);
            if isempty(isit)
                lonlatseg=[lonlatseg;lon1 lat1 lon2 lat2];  %#ok<*AGROW>
                ndseg=[ndseg;nn(j) nn(j+1)]; %#ok<*AGROW>
            end
        end
    end
end

