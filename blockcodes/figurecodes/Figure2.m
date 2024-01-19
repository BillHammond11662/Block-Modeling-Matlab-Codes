
figure(2);
clf
hold on;
PlotStates(0,0);

PlotBlocks(nodes,blocks,faults,0,0,0);
PlotNodes(nodes,blocks);
 
axis(bounds);
box on;

title('Node Numbers');