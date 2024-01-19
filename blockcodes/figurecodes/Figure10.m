
figure(10);
clf;
PlotStates(0,0);
hold on
PlotBlocks(nodes,blocks,faults,0,0,0);
axis(bounds);
xlabel('Longitude');
ylabel('Latitude');

title('Red=Data, Green=Model, Black=Residual');

s=50;

plot(lon,lat,'k.');
hv=quiver(lon(ilist),lat(ilist),s*rese/1000,s*resn/1000,0,'r');
set(hv,'LineWidth',2,'MaxHeadSize',1)

text(lon(ilist),lat(ilist),sta(ilist))
dy = diff(get(gca,'ylim'));
dx = diff(get(gca,'xlim'));

lonleg = bounds(1)+dx/8;
latleg = bounds(3)+dy/8;

hq=quiver(lonleg,latleg,s*1/1000,0,0,'k');
set(hq,'LineWidth',2,'MaxHeadSize',1)
text(lonleg,latleg+.05,'1 mm/yr');

