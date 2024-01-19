function [blocks,faults]=ReduceBlockByScore(blocks,faults,nodes,lonlatseg,score,optverb)
% [blocks,faults]=ReduceBlockByScore(blocks,faults,nodes,lonlatseg,score,optverb)
%
% Reduce the model by removing segments with score=0

if ~isempty(find(score==0, 1))
    still=1;
else
    still=0;
    disp('Not removing any by score. None with score=0');
end

if size(faults,1)==0
    blocks=[];
    return;
end

[A0,~,~,~,~,~,~,~] = blockgeoms(blocks,nodes);
Asum0=sum(A0);

while still
    
    if optverb
       disp(['Still have ' num2str(length(find(score==0))) ' to go.']);
    end
    k=find(score==0, 1);
    
    P=size(faults,1);
        
    lonlatsegf=[];
    for q=1:P
        n1=faults(q,1);
        n2=faults(q,2);
        lonlatsegf = [lonlatsegf; nodes(n1,1) nodes(n1,2)  nodes(n2,1) nodes(n2,2)];
    end

    lat1=lonlatseg(k,2);
    lon1=lonlatseg(k,1);
    lat2=lonlatseg(k,4);
    lon2=lonlatseg(k,3);
    
%     if size(lonlatsegf,2)<4
%         keyboard;
%     end
    
    %get the fault segment number 
    j=find((lonlatsegf(:,1)==lon1 & lonlatsegf(:,2)==lat1 & lonlatsegf(:,3)==lon2 & lonlatsegf(:,4)==lat2) | ...
           (lonlatsegf(:,1)==lon2 & lonlatsegf(:,2)==lat2 & lonlatsegf(:,3)==lon1 & lonlatsegf(:,4)==lat1));
    
    % if segment not present in fault list, add it
    na=[];
    nb=[];
    if isempty(j)
       na = find(nodes(:,1)==lon1 & nodes(:,2)==lat1);
       nb = find(nodes(:,1)==lon2 & nodes(:,2)==lat2);
       if isempty(na) || isempty(nb)
           disp('Uh, oh. No node numbers found for this segment');
           keyboard;
       end
       faults = [faults;[na nb nan nan]]; 
       j=size(faults,1);
    end
    
    s=Fault2Blocks(nodes,blocks,faults,j);

    % now remove fault segment j from block model
    segnds = [faults(j,1);faults(j,2)];
    [blocks2,faults2,~,~]=RemoveBlock(blocks,nodes,faults,segnds);    
    
    [A,~,~,~,~,~,~,~] = blockgeoms(blocks2,nodes);
    Asum=sum(A);
    if ((Asum0-Asum)>.001 || length(s)~=2)

        disp(['Warning: Not removing segment ' num2str(j) '.  Block removal error. Continuing.']);
        score(k)=Inf;
        faults(end,:)=[];

    else
        blocks=blocks2;
        faults=faults2;
        
        score(k)=[];
        lonlatseg(k,:)=[];
    end
    
    if isempty(find(score==0, 1))
        still=0;
    end
    
end


%%  perform imporant checks on model

CheckFaults(nodes,blocks,faults);

blocks=CheckBlocks(blocks);



