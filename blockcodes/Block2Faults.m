function flist=Block2Faults(blocks,faults,bname)
%
% flist=Block2Faults(blocks,faults,bnum)

[P,four]=size(faults);
if P<1 || four~=4
    disp('Warning: faults must be Px4');
    return;
end
    
jb=blocks.(char(bname));

if isempty(jb)
    disp(['Warning: No nodes in ' bname]);
    return;
end

% generate short fault list
flist0=[];
for i=1:length(jb)
    n=jb(i);
   ifn=find(faults(:,1)==n | faults(:,2)==n);
   flist0=[flist0; ifn];
end
flist0=unique(flist0)';

flist=[];
for k=flist0
    if length(intersect(faults(k,1:2),jb))==2
        flist = [flist;k];
    end
end
flist=unique(flist);
