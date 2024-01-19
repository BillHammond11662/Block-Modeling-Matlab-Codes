function PlotBMResultFrag(nodes,blocks,faults,nu,mu,L,omega,slip,strain,sc,cmin,cmax,opt1,opt2,opt3)
%
% PlotBMResultFrag(nodes,blocks,faults,nu,mu,L,omega,slip,strain,sc,cmin,cmax,opt1,opt2,opt3)
%
% sc is the exaggeration factor
% sc2 is rotation color scale factor (its maximum rotation rate in scale)
%
% opt1 = 1 makes patches opaque
% opt2 = = provides color coded vertical axis rotations
% opt3 = 1 is fast (no plotting of strain in blocks)

blocknames = fieldnames(blocks);
M = length(blocknames);
R = 6371009.3;

hold on;
box on;
set(gca,'dataaspectratio',[1 cosd(mean(nodes(:,2))) 1]);

if opt2==1
%     cmax = 1;
%     cmin = -1;
    cmap=colormap(jet);
    Q=length(cmap);
    m = (Q-1)/(cmax-cmin);
    omegav=rad2deg(SpinRate(blocks,nodes,omega))*1e6;  % spin rate in deg/My
end

for i=1:M
    eval(['jj=blocks.' char(blocknames(i)) ';']);

%     [vltn,vlte]=PolePred(omega(i,:)',nodes(jj,2),nodes(jj,1));
    
    % make a temporary model with just one block
    jlist = setxor(1:M,i);
    blockstemp=blocks;
    for j=1:length(jlist)
       blockstemp=rmfield(blockstemp,blocknames(jlist(j)));
    end
    
    if opt3~=1
    
        [vepred,vnpred,vrote,vrotn,vstre,vstrn,vcose,vcosn]=...
          vpredBM(faults,nodes,blockstemp,L,nu,mu,omega(i,:),slip,strain(i,:),nodes(jj,1),nodes(jj,2));
    else
      
       [vrotn,vrote]=PolePred(omega(i,:)',nodes(jj,2),nodes(jj,1));
       vstrn=zeros(size(vrotn));
       vstre=zeros(size(vrotn));
    end
    
    dlat = km2deg((vrotn+vstrn)'*sc/1000);
    dlon = km2deg((vrote+vstre)'*sc/1000);

    if opt2==1
        
        m = (Q-1)/(cmax-cmin);
        b = 1 - cmin*(Q-1)/(cmax-cmin);
        iclr=m*omegav(i) + b;
        if iclr<1
            iclr=1;
        end
        if iclr>Q
            iclr=Q;
        end
        iclr = fix(iclr);
        clr =  cmap(iclr,:);
        caxis([cmin cmax]);
        colorbar('horiz')
       
    else
        clr = [.9 .9 .9];
    end

    if opt1==1
        patch(nodes(jj,1)+dlon',nodes(jj,2)+dlat',clr);
    else
        plot(nodes(jj,1)+dlon',nodes(jj,2)+dlat','k-')
    end

end