
figure(5);
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
quiver(lon,lat,s*ve/1000,s*vn/1000,0,'r');

dy = diff(get(gca,'ylim'));
dx = diff(get(gca,'xlim'));

lonleg = bounds(1)+dx/8;
latleg = bounds(3)+dy/8;

quiver(lonleg,latleg,s*5/1000,0,0,'k');
text(lonleg,latleg+.1,'5 mm/yr');

