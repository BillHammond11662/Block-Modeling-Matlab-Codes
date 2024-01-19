function [hb,hp]=PlotBlocks(nodes,blocks,faults,opt1,opt2,opt3)
%
% hb=PlotBlocks(nodes,blocks,faults,opt1,opt2,opt3)
%
% opt1==1 draws filled patches
% opt2==1 labels blocks with names
% opt3==1 labels nodes with red node number
%
% hb are the block outline graphic objects
% hp are the block patch graphic objects

blocknames = fieldnames(blocks);
M = length(blocknames);

% PlotStates;
hold on;

set(gca,'dataaspectratio',[1 cosd(mean(nodes(:,2))) 1]);
% plon = (nodes(faults(:,1),1) + nodes(faults(:,2),1))/2;
% plat = (nodes(faults(:,1),2) + nodes(faults(:,2),2))/2;
lab=num2str((1:length(faults))');

hb=cell(M,1);
hp=cell(M,1);

for i=1:M
%     eval(['jj=blocks.' char(blocknames(i)) ';']);
    jj=blocks.(blocknames{i});
    if max(jj)>length(nodes)
        error(['node number ' num2str(max(jj)) ' is greater than length of nodes in block ' num2str(i)]);
    end
    hb{i}=plot(nodes(jj,1),nodes(jj,2),'m-','linewidth',1);
end

if opt1==1
    for i=1:M
%         eval(['jj=blocks.' char(blocknames(i)) ';']);
        jj=blocks.(blocknames{i});
        hp{i}=patch(nodes(jj,1),nodes(jj,2),[.9 .9 .9]);
    end
end

if opt2==1
    for i=1:M
%         eval(['jj=blocks.' char(blocknames(i)) ';']);
        jj=blocks.(blocknames{i});

%         if nodes(jj(1),1)==nodes(jj(end),1) && nodes(jj(1),2)==nodes(jj(end),2)
%             jj(end)=[];
%         end
        [xmn,ymn]=centroid(polyshape(nodes(jj,1),nodes(jj,2),'KeepCollinearPoints',true));
%         [GEOM,~,~]=polygeom(nodes(jj,1),nodes(jj,2));
%         xmn = GEOM(2);
%         ymn = GEOM(3);
        
        ht=text(xmn,ymn,[num2str(i) ':' char(blocknames(i))]);
        
        set(ht,'FontSize',10);
        set(ht,'HorizontalAlignment','center');
    end
end

if opt3==1
   Q=size(nodes,1);
   for i=1:Q
       text(nodes(i,1)+.01,nodes(i,2),num2str(i),'color','r');
   end
end