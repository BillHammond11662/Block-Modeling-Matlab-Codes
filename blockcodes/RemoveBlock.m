function [blocksnew,faultsnew,blocksout,faultsout]=RemoveBlock(blocks,nodes,faults,segnds)
%  [blocksnew,faultsnew,blocksout,faultsout]=RemoveBlock(blocks,nodes,faults,segnds);
%
%  given segment defined by two nodes in "segnds" this routine 
%  finds and deletes the two blocks that this segment bounds
%  and adds a new block with the combined area. 
%
%  All faults on the boundary between those blocks are removed

blocknames=fieldnames(blocks);

j=seg2blocks(blocks,segnds);
blocksout=j;

if length(j)~=2
    disp(['Warning: This segment [' num2str(segnds(1)) ' ' num2str(segnds(2)) '] is not bounded by two blocks. Exiting. ']);
    blocksnew=blocks;
    faultsnew=faults;
    faultsout=[];
    return;
end

%% check for funky block situation BEFORE

if length(j)==2
    if j(1)==j(2)
        disp('Warning: This segment funky block situation. Same block on both sides of this segment. Cancelling block removal. ');
        blocksnew=blocks;
        faultsnew=faults;
        faultsout=[];
        blocksout=[];
        return;
    end
end


%%  find node list of combined block

% keyboard;

nds1 = blocks.(blocknames{j(1)});
nds2 = blocks.(blocknames{j(2)});

x1=nodes(nds1,1);
y1=nodes(nds1,2);
x2=nodes(nds2,1);
y2=nodes(nds2,2);
[x1,y1]=poly2cw(x1,y1);
[x2,y2]=poly2cw(x2,y2);

[xu,yu] = polybool('union',x1,y1,x2,y2);
[xu, yu] = poly2ccw(xu,yu);

if any(isnan(xu)) || any(isnan(yu))
    disp('Bummer. NaNs found in xu,yu of combined block. Cancelling this block removal');
    blocksnew=blocks;
    faultsnew=faults;
    faultsout=[];
    return;
end
% figure(100);
% clf;
% plot(xu,yu,'go-','markersize',6,'linewidth',2);
% hold on;
% plot(x1,y1,'bo--','linewidth',1);
% plot(x2,y2,'r*--','linewidth',1);

ndsu = [];
for i=1:length(xu)
    %k=find(xu(i)==nodes(:,1) & yu(i)==nodes(:,2));
    % need to account for machine accuracy issue 
    k=find(abs(xu(i)-nodes(:,1))<1e-7 & abs(yu(i)-nodes(:,2))<1e-7);
    if length(k)>1
        disp('Warning: Nodes very close to one another');
    end
    if isempty(k)
       disp('doh! cant find new block node in node list.');
       keyboard;
    else
       ndsu=[ndsu;k(1)];
    end
end

blocksnew = blocks;
blocksnew=rmfield(blocksnew,blocknames{j(1)});
blocksnew=rmfield(blocksnew,blocknames{j(2)});

bnamenew = ['CB' blocknames{j(1)} '_' blocknames{j(2)}];
% matlab has a limit of 63 characters for name of field inside
% structured array
if length(bnamenew)>63
    bnamenew=bnamenew(1:63);
    if strcmp(bnamenew(end),'_')
        bnamenew(end)=[];
    end
    bnames=fieldnames(blocks);
    still = 1;
    while still
        isitunique = 1;
        for i=1:length(bnames)
            if strcmp(bnamenew,bnames{i})
                isitunique = 0;
                clast = double(bnamenew(end));
                bnamenew = [bnamenew(1:(end-1)) char(clast+1)];
            end
        end
        if isitunique==1
            still=0;
        end
    end

end
blocksnew.(bnamenew)=ndsu;

%% check to see if removal resulted in bad blocks

bnn=fieldnames(blocksnew);
j=0;
for b = bnn'
    j=j+1;
    nn=blocksnew.(char(b));
    if length(nn)<3
        disp(['Warning:  Block ' num2str(j) ': ' char(b) ' has less than 3 nodes. Not removing block.']);
        blocksnew=[];
        faultsnew=[];
        faultsout=[];
        return;
    end  
end



%% Remove faults associated with the boundary between the removed blocks

% go through all faults 
% find orphans
% remove them

P=size(faults,1);

faultsout=[];
optverb=0;
for i=1:P
   
    nf1 = faults(i,1);
    nf2 = faults(i,2);
    jf=seg2blocks(blocksnew,[nf1;nf2]);
    
    if length(jf)~=2
        if optverb
           disp(['Removing fault ' num2str(i)]);
        end
        faultsout=[faultsout;i];
    end
    
end
faultsnew=faults;
faultsnew(faultsout,:)=[];

%% check for funky block situation AFTER

ndschk = blocksnew.(bnamenew);

% truncate end node if its there
if ndschk(1)==ndschk(end)
    ndschk(end)=[];
end

for i=1:length(ndschk)
    ndc1 = ndschk(i);
    if i==length(ndschk)
        ndc2 = ndschk(1);
    else
        ndc2 = ndschk(i+1);
    end
    
    btest=seg2blocks(blocks,[ndc1;ndc2]);
    if length(btest)==2
        if btest(1)==btest(2)
            disp('Warning! Cancelling this block combination since it results in a funky block!');
            blocksnew=blocks;
            faultsnew=faults;
            faultsout=[];
            blocksout=[];
        end
    end
end

