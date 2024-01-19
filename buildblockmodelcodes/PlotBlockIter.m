
clf;

A=gcf;
if (optverb==1)
    disp(['Plotting on Figure ' num2str(A.Number)]);
end

PlotBlocks(nodes,blocks,faults,1,1,1);
hold on
PlotFaults(blocks,nodes,faults,0,0,1);
set(gca,'dataaspectratio',[1 cosd(median(nodes(:,2))) 1])
axis(bnds);
title([num2str(length(fieldnames(blocks))) ' Blocks and ' num2str(size(faults,1)) ' Faults']);
