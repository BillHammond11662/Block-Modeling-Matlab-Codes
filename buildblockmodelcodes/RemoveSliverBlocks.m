function [blocks,faults]=RemoveSliverBlocks(blocks,faults,nodes,bnds,degsmall,optverb)
% [blocks,faults]=RemoveSliverBlocks(blocks,faults,nodes,bnds,degsmall,optverb)
%
% "Sliver" blocks are those that have very small interior angles and also
% lie along the boundary. Because they contribute nothing to the slip rate
% estimation, they are removed completely rather than combining with
% adjacent blocks.

[~,~,~,~,~,~,amin,~] = blockgeoms(blocks,nodes);

bnames = fieldnames(blocks);
bnamesout = bnames(amin<degsmall);

for b=bnamesout'
    
    if isfield(blocks,char(b))
        
        nn = blocks.(char(b));
        if nn(1)==nn(end)
            nn(end)=[];
        end
        
        angs  = blockangles(nodes(nn,1),nodes(nn,2));
        
        isout=false;
        for j=1:length(nn)
            lon1=nodes(nn(j),1);
            lat1=nodes(nn(j),2);
            if j==length(nn)
                lon2=nodes(nn(1),1);
                lat2=nodes(nn(1),2);
            else
                lon2=nodes(nn(j+1),1);
                lat2=nodes(nn(j+1),2);
            end
            
            if ((lon1==bnds(1) && lon2==bnds(1)) || ...
                    (lon1==bnds(2) && lon2==bnds(2)) || ...
                    (lat1==bnds(3) && lat2==bnds(3)) || ...
                    (lat1==bnds(4) && lat2==bnds(4)) )
                
                if j<length(nn)
                    if (angs(j)<degsmall || angle(j+1)<degsmall)
                        isout=true;
                    end
                else
                    if (angs(j)<degsmall || angle(1)<degsmall)
                        isout=true;
                    end
                end
                
            end
            
        end
        
        if isout
            if optverb
               disp(['Removing block: ' char(b)]);
            end
            fnums=Block2Faults(blocks,faults,b);
            blocks=rmfield(blocks,b);
            faults(fnums,:)=[];
        end
        
    end
end

%%  perform imporant checks on model

CheckFaults(nodes,blocks,faults);

blocks=CheckBlocks(blocks);
