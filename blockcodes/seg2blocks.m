function blist=seg2blocks(blocks,segnds)
%
% blist=seg2blocks(blocks,segnds)
%
% finds blocks associated with a pair of nodes.
% nodes must be adjacent to specify valid segment.

% check to make sure segnds is 2x1
[n,m]=size(segnds);
if ~(n==2 && m==1)
    error('segnds must be 2x1');
end

% write something here to check to make sure nodes are adjacent
%
   
blocknames=fieldnames(blocks);
M=length(blocknames);

blist = [];
for i=1:M
    nds=blocks.(blocknames{i});
    
    % remove wrap-around node if there is one
    if nds(1)==nds(end)
        nds(end)=[];
    end
    
    for j=1:(length(nds)-1)
        if ((segnds(1)==nds(j) && segnds(2)==nds(j+1)) || (segnds(2)==nds(j) && segnds(1)==nds(j+1)))
            blist=[blist;i]; %#ok<AGROW>
        end
    end
    if ((segnds(1)==nds(1) && segnds(2)==nds(end)) || (segnds(2)==nds(1) && segnds(1)==nds(end)))
        blist=[blist;i]; %#ok<AGROW>
    end
end