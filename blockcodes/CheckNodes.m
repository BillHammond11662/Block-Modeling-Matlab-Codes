function [flag,isus]=CheckNodes(nodes,faults,blocks)

disp('Checking for Duplicate Nodes')

flag=0;

eps = km2deg(.1);
blnames = fieldnames(blocks);
M=length(blnames);

isus=[]; % list of suspect nodes

N=size(nodes,1);
for i=2:N
    imnot = 1:N;
    imnot(imnot>=i)=[];
    [D,~]=baz(nodes(i,2),nodes(i,1),nodes(imnot,2),nodes(imnot,1));
    [minD,imin]=min(D);
    f1=0;
    if minD==0
        flag=1;
        disp(['Node ' num2str(i) ' is same as node(s) ' num2str(imnot(imin))]);
        f1=1;
    elseif minD<eps
        flag=1;
        disp(['Node ' num2str(i) ' is extremely close to node(s) ' num2str(imnot(imin))]);
        f1=1;
    end

    % does it matter? is this node in faults/blocks?
    if f1==1
        usedinfaults=1;
        if isempty(find(faults(:,1:2)==i, 1))
            usedinfaults=0;
        end
        usedinblocks=0;
        for j=1:M
            nlist=blocks.(blnames{j});
            if ~isempty(intersect(nlist,i))
                usedinblocks=1;
            end
        end
        if (~usedinfaults && ~usedinblocks)
            disp('   But it doesn''t matter because its not used in faults of blocks.')
            f1=0;
        else
            isus=[isus;i;imnot(imin)];
        end
    end
end

if flag==1
    disp('Found duplicate nodes.');
end