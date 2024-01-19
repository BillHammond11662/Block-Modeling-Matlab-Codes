function j=Fault2Blocks(nodes,blocks,faults,k)
% j=Fault2Blocks(nodes,blocks,faults,k);
%
% nodes,blocks,faults are standard block model
% k fault number
%
% output j is indeces of
%  the two blocks adjacent to a given fault


blocknames = fieldnames(blocks);
M = length(blocknames);

fn1=faults(k,1);
fn2=faults(k,2);

% shorten the list of potential blocks for speed since the "intersect"
% command below is time consuming.
blist=[];
for i=1:M
    nds=Block2Nodes(blocks,i);
    if ~isempty(find(nds==fn1,1)) || ~isempty(find(nds==fn2,1))
        blist=[blist;i];
    end
end

j=[];
for jj=blist'
    jb=Block2Nodes(blocks,jj);
    if (length(intersect([fn1 fn2],jb))==2)
        j=[j;jj];
    end
end

j2=[];
tiny=.0001;
if length(j)>2

    % find fault midpoint
    px = (nodes(fn1,1) + nodes(fn2,1))/2;
    py = (nodes(fn1,2) + nodes(fn2,2))/2;

    % find elements (blocks) of j that have midpoint

    for kk=1:length(j)
        %         eval(['jb=blocks.' char(blocknames(j(kk))) ';']);
        jb=Block2Nodes(blocks,j(kk));
        in1=inpolygon(px+tiny,py+tiny,nodes(jb,1),nodes(jb,2));
        in2=inpolygon(px-tiny,py-tiny,nodes(jb,1),nodes(jb,2));

        if (in1 || in2)
            j2 = [j2;j(kk)];
        end
    end

    j = j2;

end

if length(j)~=2
    %disp(['Warning: Could not place fault ' num2str(k) ' length of j=' num2str(length(j2))]);
else
    
    % decide which block is on 'left' of fault
    % if j(1) is on right switch the order so we always subtract
    % the block on the right from the one on the left.
    
    % eval(['jj=blocks.' char(blocknames(j(1))) ';']);
    
    jj=Block2Nodes(blocks,j(1));
    inn = faults(k,1:2);
    mtch=0;
    for q = 1:(length(jj)-1)
        if all(inn==jj(q:q+1))
            mtch=1;
        end
    end
    if mtch==0
        j = flipud(j);
    end
end
