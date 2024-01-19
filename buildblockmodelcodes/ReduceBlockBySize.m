function [blocks,nodes,faults]=ReduceBlockBySize(blocks,nodes,faults,Amin,bnds,optverb)
% [blocks,nodes,faults]=ReduceBlockBySize(blocks,nodes,faults,Amin,bnds,optverb)
% 
% Remove blocks in model that are smaller than Amin.
%
% Routine combines small block with adjacent block touching longest edge of small block  
% that is not a fault or model boundary. 

%  Identify small blocks 

[A,~,~,~,~,~,~,~] = blockgeoms(blocks,nodes);

% remove each block by taking its longest side that is not a fault segment
bnames = fieldnames(blocks);
bnamesout = bnames(A<Amin);

if isempty(bnamesout)
    disp(' ');
    disp(['All blocks have area greater than ' sprintf('%.2f',Amin) '. Nothing to do.']);
    return;
end

for b=bnamesout'
    
    if optverb
       disp(['Trying to remove small block #: ' char(b)]);
    end
    
    if isfield(blocks,char(b))
        nn =eval(['blocks.' char(b)]);
        if nn(1)==nn(end)
            nn(end)=[];
        end
        
        % find candidate segments to remove
        isfs=false(size(nn));
        L=nan(size(nn));
        for j=1:length(nn)
            ns1=nn(j);
            lon1=nodes(ns1,1);
            lat1=nodes(ns1,2);
            if j==length(nn)
                ns2=nn(1);
                lon2=nodes(ns2,1);
                lat2=nodes(ns2,2);
            else
                ns2=nn(j+1);
                lon2=nodes(ns2,1);
                lat2=nodes(ns2,2);
            end
            [L(j),~] = distance(lat1,lon1,lat2,lon2);
            
            % is this segment a boundary segment?
            % there is a problem here in that previous iterations could
            % remove blocks with an edge on the boundary.  Thus this condition might not work to identify a bondary block
            % replace this with "seg2blocks()"
%             if ((lon1==bnds(1) && lon2==bnds(1)) || ...
%                     (lon1==bnds(2) && lon2==bnds(2)) || ...
%                     (lat1==bnds(3) && lat2==bnds(3)) || ...
%                     (lat1==bnds(4) && lat2==bnds(4)) )
%                 isfs(j)=true;
%             end
            js = seg2blocks(blocks,[ns1;ns2]);
            if length(js)<2
                isfs(j)=true;
            end
            
            % is this segment a fault segment?
            fnum=[];
            if j==length(nn)
                fnum=nodes2fault(faults,nn(j),nn(1));
            else
                fnum=nodes2fault(faults,nn(j),nn(j+1));
            end
            if ~isempty(fnum)
                isfs(j)=true;
            end
            
        end
        
        if all(isfs)
            
            disp(['No candidate segments. Not removing small block #: ' char(b)]);
            
        else
            
            % find the longest segment with isfs false
            [~,isrt]=sort(L,'descend');
            still=1;
            k=1;
            isegrm=nan;
            while still
                if ~isfs(isrt(k))
                    isegrm=isrt(k);
                    still=0;
                end
                k=k+1;
            end
            
            if isnan(isegrm)
                if optverb
                    disp(['   Not removing small block #: ' char(b)]);
                end
            else
                if optverb
                    disp(['   Actually removing small block #: ' char(b)]);
                end
                
                if isegrm==length(nn)
                    segnds=[nn(isegrm);nn(1)];
                else
                    segnds=[nn(isegrm);nn(isegrm+1)];
                end
                
                [blocks,faults,~]=RemoveBlock(blocks,nodes,faults,segnds);
                
            end
            
        end
    else
        if optverb
           disp(['Cannot find block ' char(b) ', may have been combined already']);
        end
    end
end

%%  perform imporant checks on model

CheckFaults(nodes,blocks,faults);

blocks=CheckBlocks(blocks);



