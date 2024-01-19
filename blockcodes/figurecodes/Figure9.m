
figure(9);
clf;

PlotStates(0,0);
hold on;

if isempty(vepred)
    return;
end

% If block is blank means there was not enough GPS stations on block to get
% an RMS misfit
[chi2dof,dof,n,rmsn,rmse,rmstot]=...
    BlockMisfit(blocks,nodes,lon(ilist),lat(ilist),ve(ilist),vn(ilist),se(ilist),sn(ilist),vepred(ilist)*1000,vnpred(ilist)*1000,0);

cmin = 0;
cmax = 3;
cmap=colormap(jet);
Q=length(cmap);
m = (Q-1)/(cmax-cmin);
b = 1 - cmin*(Q-1)/(cmax-cmin);

M=length(fieldnames(blocks));

for i=1:M
%     eval(['jj=blocks.' char(blocknames(i)) ';']);
    jj=blocks.(blocknames{i});
    if ~isnan(rmstot(i))
        iclr=m*rmstot(i) + b;
        if iclr<1
            iclr=1;
        end
        if iclr>Q
            iclr=Q;
        end
        iclr = fix(iclr);
        clr =  cmap(iclr,:);
        
        patch(nodes(jj,1),nodes(jj,2),clr);
    end
    
end

clim([cmin cmax]);
colorbar('horiz')
box on;
axis(bounds);
hold on;
title('RMS Misfit per block');
set(gca,'dataaspectratio',[1 cosd(40) 1]);

