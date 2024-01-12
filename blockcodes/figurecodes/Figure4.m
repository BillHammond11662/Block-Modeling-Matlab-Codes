

figure(4);
clf
hold on;
PlotBlocks(nodes,blocks,faults,0,1,0);
hp=plot(lon,lat,'g.');
set(hp,'color',[0 .5 0]);
set(hp,'MarkerSize',10)
axis(bounds);
box on;

