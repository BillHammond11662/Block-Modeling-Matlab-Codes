
figure(7);
clf;

PlotBlocks(nodes,blocks,faults,0,0,0);
hold on;
sc=3e6;
cmin=-1;
cmax=1;
PlotBMResultFrag(nodes,blocks,faults,nu,mu,L,omega,slip,strain,sc,cmin,cmax,1,1,0)
axis(bounds + [-.2 0 0 .2]);

title({'Rotation/Vertical Axis Spin Rate',...
      ['RMS East=' sprintf('%.2f',rmse) ' mm/yr,   RMS North=' sprintf('%.2f',rmsn) ' mm/yr'],...
       'Blue is clockwise, Red is counter-clockwise'});
