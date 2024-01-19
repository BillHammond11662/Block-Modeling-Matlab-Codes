
%figure(pnum);
clf;

PlotFaults(blocks,nodes,faults,0,0,1);
hold on

[G,~]=size(lonlatseg);

for i=1:G
    
    lon1=lonlatseg(i,1);lat1=lonlatseg(i,2);
    lon2=lonlatseg(i,3);lat2=lonlatseg(i,4);
    
    if score(i)==1
        plot([lon1 lon2],[lat1 lat2],'m-','linewidth',2);
    elseif score(i)==2
        plot([lon1 lon2],[lat1 lat2],'r-','linewidth',2);
    elseif score(i)==10
        plot([lon1 lon2],[lat1 lat2],'g-','linewidth',2);
    elseif score(i)==Inf
        plot([lon1 lon2],[lat1 lat2],'b-','linewidth',2);
    else
        plot([lon1 lon2],[lat1 lat2],'-','linewidth',1,'color',[.5 .5 .5]);
    end
        
end
axis(bnds);
set(gca,'dataaspectratio',[1 cosd(38) 1])

title({[num2str(length(fieldnames(blocks))) ' Blocks and ' num2str(size(faults,1)) ' Faults'],...
    'Gray=0 - Magenta=1 - Red=2 - Green=10 - Blue=Inf'});