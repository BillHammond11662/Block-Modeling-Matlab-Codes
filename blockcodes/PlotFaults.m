function PlotFaults(blocks,nodes,faults,opt,opt2,opt3)
%
%  PlotFaults(blocks,nodes,faults,opt,opt2,opt3);
%
% if opt  = 1 plot nodes too
% if opt2 = 1 plot dip markers
% if opt3 = 1 plot fault numbers

%PlotStates(0,0);
if isempty(faults)
    return;
end
hold on;

% plot blocks
if ~isempty(blocks)
    blocknames = fieldnames(blocks);
    M = length(blocknames);
    for i=1:M
        eval(['jj=blocks.' char(blocknames(i)) ';']);
        plot(nodes(jj,1),nodes(jj,2),'m-');
    end
else
    M=0;
end

% plot nodes
if opt==1
    plot(nodes(:,1),nodes(:,2),'k.');
end

% plot faults
plon = (nodes(faults(:,1),1) + nodes(faults(:,2),1))/2;
plat = (nodes(faults(:,1),2) + nodes(faults(:,2),2))/2;
lab=num2str((1:length(faults))');
sc = .1;
for i=1:size(faults,1)
    plot([nodes(faults(i,1),1) nodes(faults(i,2),1)],...
        [nodes(faults(i,1),2) nodes(faults(i,2),2)],'g-','linewidth',3);
    if opt3
       text(plon(i),plat(i),lab(i,:))
    end


    xmid = (nodes(faults(i,1),1) + nodes(faults(i,2),1))/2;
    ymid = (nodes(faults(i,1),2) + nodes(faults(i,2),2))/2;


    dip = faults(i,4);

    dx = nodes(faults(i,2),1) - nodes(faults(i,1),1);
    dy = nodes(faults(i,2),2) - nodes(faults(i,1),2);

    % rotate 90? CCW to point towards positive dip direction
    % and normalize then scale by sc
    L = norm([dx;dy]);

    if opt2

        dxp= -sc*dy*cosd(dip)*sign(dip)/L;
        dyp=  sc*dx*cosd(dip)*sign(dip)/L;

        plot([xmid xmid+dxp],[ymid ymid+dyp],'k-','linewidth',3);

    end
end

set(gca,'dataaspectratio',[1 cosd(mean(nodes(:,2))) 1]);