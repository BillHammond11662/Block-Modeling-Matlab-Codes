

Fn=gcf;
figure(Fn.Number+1);
clf
PlotStates(0,0);
hold on;
for i=1:Nf
    plot(geol(i).lonseg,geol(i).latseg,'r-','linewidth',3);
end
[hb,hp]=PlotBlocks(nodes,blocks,faults,0,1,0);
for i=1:length(hb)
    handle = hb{i};
    handle.Color=[0 .5 0];
end
axis(bounds);
xlabel('Longitude');
ylabel('Latitude');
title('Block Boundaries (green), Input Faults (red), and Names')

Fn=gcf;
figure(Fn.Number+1);
clf
hold on;
PlotStates(0,0);
PlotBlocks(nodes,blocks,faults,0,0,0);
PlotNodes(nodes,blocks);
axis(bounds);
box on;
title('Node Numbers');


Fn=gcf;
figure(Fn.Number+1);
clf;
title('Fault Numbers')
PlotStates(0,0);
hold on;
PlotFaults(blocks,nodes,faults,1,1,1);
title('Fault Numbers and Dip Markers');
axis(bounds);