function blocks = CheckBlocks(blocks)
% blocks = CheckBlocks(blocks)
% 

bnames=fieldnames(blocks);
j=0;
for b = bnames'
    j=j+1;
    %eval(['nn=blocks.' char(b) ';']);
    nn=blocks.(char(b));
    in=find(nn(1)==nn);
    if length(in)>2
        disp(['Warning: More than 2 instances of first node in block ' char(b)]);
        disp('Trying to fix it');
        if nn(end)==nn(end-1)
            eval(['blocks.' char(b) ' = nn(1:end-1);']);
            disp('Succeded');
        end
        if nn(1)==nn(2)
            eval(['blocks.' char(b) ' = nn(2:end);']);
            disp('Succeded');
        end
    end
    
    %
    if length(nn)<3
        disp(['Warning:  Block ' num2str(j) ': ' char(b) ' has less than 3 nodes.']);
        %         keyboard;
        
    end
    
%     if nn(1)==nn(end)
%         nn2=nn(1:end-1);
%     else
%         nn2=nn;
%     end
%     polygeom(nodes(nn2,1),nodes(nn2,2));

end


%% add section to check for funky block.  See RemoveBlock2 for code


