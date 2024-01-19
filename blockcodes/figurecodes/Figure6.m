
figure(6);
clf;
sigplot=1;
lenleg=3;
fflag=0;
gmtflag=0;
mode=1;
box on;

rmse=1000*rms((ve(ilist)/1000)-vepred(ilist));
rmsn=1000*rms((vn(ilist)/1000)-vnpred(ilist));

PlotStates(0,0);
hold on;

ts = 3; % adjusts width/fatness of line designating strike slip rates

PlotBMResult(nodes,blocks,faults,slip,s_slip,...
    sigplot,latleg,lonleg,lenleg,fflag,gmtflag,ts);

axis(bounds);

blocknames = fieldnames(blocks);

if any(any(~isnan(strain)))
    s=1;
    for j=1:size(strain,1)
        
        if any(abs(strain(j,:))>s_strain(j,:))
            
            e_phiphi = strain(j,1);
            e_thetaphi = strain(j,2);
            e_thetatheta = strain(j,3);
            s_e_phiphi = s_strain(j,1);
            s_e_thetaphi = s_strain(j,2);
            s_e_thetatheta = s_strain(j,3);
            nds = blocks.(char(blocknames(j)));
            
            [geom,~,~] = polygeom( nodes(nds,1), nodes(nds,2) );
            lon0 = geom(2);
            lat0 = geom(3);
            
            [alphaaz,e1,e2,se1,se2]=plot_tensor(lat0,lon0,e_phiphi,e_thetatheta,...
                e_thetaphi,s_e_phiphi,s_e_thetatheta,s_e_thetaphi,s,1,1,3,'r','b');
            
        end
        
    end
    
end

title(['Slip Rate. # Blocks=' num2str(length(blocknames)) ', RMS East=' num2str(rmse) ',   RMS North=' num2str(rmsn)]);

