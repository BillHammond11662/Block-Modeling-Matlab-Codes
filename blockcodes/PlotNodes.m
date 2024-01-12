function PlotNodes(nodes,blocks)
%
% PlotNodes(nodes,blocks);
%

blocknames = fieldnames(blocks);
M = length(blocknames);

hold on;

set(gca,'dataaspectratio',[1 cosd(mean(nodes(:,2))) 1]);

for i=1:M
    eval(['jj=blocks.' char(blocknames(i)) ';']);
    plot(nodes(jj,1),nodes(jj,2),'m-')
end

plot(nodes(:,1),nodes(:,2),'k.');
[N,dum]=size(nodes);

for i=1:N
   ht=text(nodes(i,1)+.01,nodes(i,2),num2str(i));
   set(ht,'color','r');
end;
