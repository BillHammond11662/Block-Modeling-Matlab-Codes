
figure(1);
clf
PlotStates(0,0);
hold on;
PlotBlocks(nodes,blocks,faults,0,1,0);
axis(bounds);
xlabel('Longitude');
ylabel('Latitude');
title('Block Boundaries and Names')

