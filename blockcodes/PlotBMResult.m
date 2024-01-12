function PlotBMResult(nodes,blocks,faults,slip,s_slip,sigplot,latleg,lonleg,lenleg,fflag,gmtflag,ts)
%
%  PlotBMResult(nodes,blocks,faults,slip,s_slip,sigplot,latleg,lonleg,lenleg,fflag,gmtflag,ts);
%
% sigplot = 0 says do NOT plot insignificant extensional slip
% sigplot = 1 says do.
%
% latleg,lonleg are coordinates of legend
% lenleg is length of bar in mm/yr
%
% fflag=1 plots the fault numbers too
% ts = is slip rate plotting scale factor
%
% gmt flag outputs files for gmtplotting: 'normal.txt' 'thrust.txt'
% 'dextral.txt' 'sinestral.txt'

if gmtflag
    fid1=fopen('normal.txt','w');
    fid2=fopen('thrust.txt','w');
    fid3=fopen('dextral.txt','w');
    fid4=fopen('sinestral.txt','w');
    if (fid1==0 || fid2==0 || fid3==0 || fid4==0)
        error('Problem opening gmt output files for writing.');
    end
end
gmtfact = 3;    

if ~isempty(slip)

    blocknames = fieldnames(blocks);
    M = length(blocknames);

    PlotStates(0,0)
    hold on;
    axis([-124 -113 36 44]);
    set(gca,'dataaspectratio',[1 cosd(40) 1]);
    plon = (nodes(faults(:,1),1) + nodes(faults(:,2),1))/2;
    plat = (nodes(faults(:,1),2) + nodes(faults(:,2),2))/2;
    lab=num2str((1:length(faults))');
    for i=1:M
        eval(['jj=blocks.' char(blocknames(i)) ';']);
        plot(nodes(jj,1),nodes(jj,2),'m-')
    end
    
    for i=1:size(faults,1)

        if slip(i,1)>0
            plot([nodes(faults(i,1),1) nodes(faults(i,2),1)],...
                [nodes(faults(i,1),2) nodes(faults(i,2),2)],'r-','linewidth',slip(i,1)*1000*ts);

            if gmtflag
                fprintf(fid4,'%9.4f %9.4f %8.4f \n',[nodes(faults(i,1),1)  nodes(faults(i,1),2)   gmtfact*abs(slip(i,1))*1000*ts ]);
                fprintf(fid4,'%9.4f %9.4f %8.4f \n>\n',[nodes(faults(i,2),1)  nodes(faults(i,2),2)   gmtfact*abs(slip(i,1))*1000*ts ]);
            end

        else
            plot([nodes(faults(i,1),1) nodes(faults(i,2),1)],...
                [nodes(faults(i,1),2) nodes(faults(i,2),2)],'k-','linewidth',-slip(i,1)*1000*ts);
            
            if gmtflag
                fprintf(fid3,'%9.4f %9.4f %8.2f \n',[nodes(faults(i,1),1)  nodes(faults(i,1),2)  gmtfact*abs(slip(i,1))*1000*ts ]);
                fprintf(fid3,'%9.4f %9.4f %8.2f \n>\n',[nodes(faults(i,2),1)  nodes(faults(i,2),2)  gmtfact*abs(slip(i,1))*1000*ts]);
            end

        end
        if fflag
            text(plon(i)+.1,plat(i),lab(i,:),'fontsize',9);
        end

        % plot the horizontal projection of bk

        if (abs(slip(i,2)/s_slip(i,2))>2 || sigplot==1)

            proj = cosd(median(nodes(:,2)));

            dy = nodes(faults(i,2),2) - nodes(faults(i,1),2);
            dx = nodes(faults(i,2),1) - nodes(faults(i,1),1);
            r = sqrt(dy.^2 + dx.^2);
            dy = dy/r;
            dx = proj*dx/r;

            scal = .025*slip(i,2)*1000*cosd(faults(i,4))*ts;

            if scal<0
%                 h=quiver(plon(i),plat(i),-dy*scal,dx*scal,0,'b');
%                 h2=quiver(plon(i),plat(i),dy*scal,-dx*scal,0,'b');
                
                plot([plon(i)-dy*scal plon(i)+dy*scal],[ plat(i)+dx*scal  plat(i)-dx*scal], 'b-','linewidth',3);
                                
                if gmtflag
                    fprintf(fid1,'%9.4f %9.4f \n',[plon(i)-dy*scal plat(i)+dx*scal]);
                    fprintf(fid1,'%9.4f %9.4f \n>\n',[[plon(i)+dy*scal plat(i)-dx*scal]]);
                end
                
            else
%                 h=quiver(plon(i),plat(i),-dy*scal,dx*scal,0,'c');
%                 h2=quiver(plon(i),plat(i),dy*scal,-dx*scal,0,'c');
                
                plot([plon(i)-dy*scal plon(i)+dy*scal],[ plat(i)+dx*scal  plat(i)-dx*scal], 'c-','linewidth',3);

                
                if gmtflag
                    fprintf(fid2,'%9.4f %9.4f \n',[plon(i)-dy*scal plat(i)+dx*scal]);
                    fprintf(fid2,'%9.4f %9.4f \n>\n',[plon(i)+dy*scal plat(i)-dx*scal]);
                end
                
            end
%             set(h,'linewidth',3);
%             set(h,'showarrowhead','off');
%             set(h2,'linewidth',3);
%             set(h2,'showarrowhead','off');
        end

    end

    if (nargin>6 && ~isempty(lonleg))
        scal = .10*lenleg;
        h=quiver(lonleg,latleg-scal/2,0,scal,0,'b');
        set(h,'linewidth',3);
        set(h,'showarrowhead','off');
        
        plot([lonleg-scal*proj lonleg+scal*proj],...
                [latleg latleg],'k-','linewidth',lenleg*ts);
        ht=text(lonleg,latleg+scal*1.1,{[num2str(lenleg) ' mm/yr'],'Normal/Dextral'});
        set(ht,'HorizontalAlignment','center');
    end

%    plot(plon,plat,'wo');

else
    disp('No Slip Information.');
end

if gmtflag
    fclose(fid1);
    fclose(fid2);
    fclose(fid3);
    fclose(fid4);
end