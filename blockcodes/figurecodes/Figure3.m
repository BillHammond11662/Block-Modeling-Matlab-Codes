

figure(3);
clf;

title('Fault Numbers')
PlotStates(0,0);
hold on;

PlotFaults(blocks,nodes,faults,1,1,1);
title('Fault Numbers and Dip Markers');

axis(bounds);