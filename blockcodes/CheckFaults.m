function [flag]=CheckFaults(nodes,blocks,faults)
% [flag]=CheckFaults(nodes,blocks,faults)
%
% Goes through each fault and checks if there are duplicates
%
%

disp(' ');
disp('Checking for duplicate faults...');

[P,four] = size(faults);
if four~=4
    error('faults is not Px4');
end

flag=0;
for i=1:P
   
    n1 = faults(i,1);
    n2 = faults(i,2);
    
    i1=find(faults(:,1)==n1 & faults(:,2)==n2);
    i2=find(faults(:,1)==n2 & faults(:,2)==n1);
    
    inot=setxor(i,[i1;i2]);
    
    if ~isempty(inot)
        disp(['Fault #' num2str(i) ' duplicate with faults: ']);
        for j=1:length(inot)
            disp(['   ' num2str(inot(j))]);
        end
        flag=1;
    end
end

if flag==0
    disp('  No duplicate faults!');
else
    disp('  There were duplicate faults.');
end

%% 

disp('Checking for very shallow dips...');

idp=find(abs(faults(:,4))<=1);
f2=0;
if ~isempty(idp)
    for i=1:length(idp)
        disp(['Fault ' num2str(idp(i)) ' has too shallow a dip.']);
        flag=1;
        f2=1;
    end
end

if f2==0
   disp('  No faults with too shallow a dip.');  
else
   disp('  There are faults with too-shallow of a dip.');
end

%% 

disp('Checking for blocks and faults with non-existent node numbers....');

maxnode = max([faults(:,1);faults(:,2)]);
[N,~]=size(nodes); 
f3=0;

if maxnode>N
    disp(' Max node number is greater than number of nodes!');
    im1=find(maxnode==faults(:,1));
    im2=find(maxnode==faults(:,2));
    im=unique([im1;im2]);
    for j=1:length(im)
       disp(['Node ' num2str(im(j)) ' in faults is not in node list.']);
       flag=1;
       f3=1;
    end
    return;
end

if f3==0
   disp('  All called nodes exist!');    
end

%% 

disp('Checking for faults not associated with blocks...');
f4=0;

for i=1:P
    
   j=Fault2Blocks(nodes,blocks,faults,i);
   
   if length(j)~=2
       disp(['   Fault ' num2str(i) ' may be illegal. It is associated with ' num2str(length(j)) ' blocks.']);
       flag=1;
       f4=1;
   end
   
end

if f4==0
   disp('  No illegal faults!');    
end

%%

disp('Checking for faults with length==0')
f5=0;
for i=1:P
    n1=faults(i,1);
    n2=faults(i,2);
    [flen,~]=baz(nodes(n1,2),nodes(n1,1),nodes(n2,2),nodes(n2,1));
    if flen==0
        disp(['  fault ' num2str(i) ' has zero length.']);
        flag = 1;
        f5=1;
    end
end

if f5==0
   disp('  No zero length faults!');    
end

%% any NaN's in the fault matrix?

if any(any(find(isnan(faults))))
   
    disp('Warning: There are NaN''s in the fault matrix.');
    flag=1;

end

%% Any faults have zero dip?

if any(find(isnan(faults(:,4)==0)))
    disp('Warning: There are faults with zero dip.')
end
