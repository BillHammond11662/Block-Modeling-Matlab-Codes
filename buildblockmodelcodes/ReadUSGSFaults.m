
% clear

% citation
% Hatem, A.E., Collett, C.M., Briggs, R.W., Gold, R.D., Angster, S.J., Powers, P.M., Field, E.H., Anderson, M., 
% Ben-Horin, J.Y., Dawson, T., DeLong, S., DuRoss, C., Thompson Jobe, J., Kleber, E., Knudsen, K.L., Koehler, R., 
% Koning, D., Lifton, Z., Madin, I., Mauch, J., Morgan, M., Pearthree, P., Pollitz, F., Scharer, K., Sherrod, B., 
% Stickney, M., Wittke, S., and Zachariasen, J., 2022, Earthquake geology inputs for the U.S. National Seismic 
% Hazard Model (NSHM) 2023 (western US) (ver. 2.0, February 2022): U.S. Geological Survey data release, 
% https://doi.org/10.5066/P9AU713N.
%
% obtained at https://www.sciencebase.gov/catalog/item/5fe3bbbcd34ea5387deb4a85

% wdir = 'wherever you have the USGS database'; 
wdir = '~/Science/Datasets/FaultDataUSGS/NSHMP2023_v2/NSHM23_FSD';

file = [wdir '/NSHM23_FSD_v2.shp'];
S=shaperead(file);

lonmin=-119;
lonmax=-117;
latmin=35.5;
latmax=37.5;
bnds = [lonmin lonmax latmin latmax];


%%

disp(' ');
disp('Plotting faults in database...');

figure(1);
clf;
PlotStates(0,0);
hold on;
plot(bnds([1 2 2 1 1]),bnds([3 3 4 4 3]),'r-');

minx=180;
maxx=-180;
miny=90;
maxy=-90;
for i=1:length(S)
    
    x=S(i).X;
    y=S(i).Y;
    plot(x,y,'b-');
    
    if max(x)>maxx
        maxx=max(x);
    end
    if min(x)<minx
        minx=min(x);
    end 
    if max(y)>maxy
        maxy=max(y);
    end
    if min(y)<miny
        miny=min(y);
    end 
    
end
axis([minx maxx miny maxy]);
set(gca,'dataaspectratio',[1 cosd(40) 1]);
drawnow;
    
%%

disp(' ');
disp('Generating geol structured array for block modeling...');

geol=S;
Nf = length(S);

for i=1:Nf
    
    % strip the NaNs off the end of the fault segs.
    lon = geol(i).X;
    lat = geol(i).Y;
    lon(isnan(lon))=[];
    lat(isnan(lat))=[];
    geol(i).lonseg = lon;
    geol(i).latseg = lat;
    geol(i).dip = geol(i).DipDeg;
    geol(i).dip_dir=geol(i).DipDir;

    klist=[];
    for k=1:length(SS)
        if str2double(S(i).FaultID) == str2double(SS(k).FaultID)
            klist=[klist;k];
        end
    end
    
    geol(i).SlipRateRecord = klist;
end

geol=rmfield(geol,'X');
geol=rmfield(geol,'Y');
geol=rmfield(geol,'DipDeg');
geol=rmfield(geol,'DipDir');


%% make edits to fault database

% exclude all faults completely outside the domain of interest
gout = [];
for i=1:Nf
    lonlatbox = geol(i).BoundingBox;
    P1=polyshape(lonlatbox([1 2 2 1 1],1),lonlatbox([1 1 2 2 1],2));
    P2=polyshape(bnds([1 2 2 1 1]),bnds([3 3 4 4 3]));
    if area(intersect(P1,P2))==0
        gout=[gout;i];
    end
end
geol(gout)=[];

% for i=1:length(geol)
%     str = [num2str(i) ': '  num2str(geol(i).FaultID) ' ' geol(i).FaultName];
%     disp(str)
% end


%%

% Remove sections of faults
% This was the painstaking part, where visual review of each fault was
% needed to makes sure sections that crossed other faults where removed.

rmseg = {'304',12;
         '2003',8;
         '157',9;
         '75',6;
         '281',3;
         '190',5;
         '145',4:5;
         '271',1;
         '266',1:2;
         '260',1:3;
         '229',1;
         '702',7;
         '168',6;
         '704',1;
         '338',1;
         '335',1;
         '224',11;
         '68',1;
         '58',1:2;
         '195',14;
         '151',9;
         '191',1;
         '186',13;
         '180',1;
         '329',1;
         '327',22:25;
         '59',1;
         '1001',7;
         '62',1;
         '2508',8};

% temporary to check and see if these removals are still necessary given
% the new version of the database.

figure(1);
hold on;

[Rf,~]=size(rmseg);
% n=nan(size(geol));
% for i=1:length(geol)
%     disp([num2str(i) ': ' num2str(geol(i).FaultID)]);
%     n(i)=geol(i).FaultID;
% end

for i=1:Rf
    for k=1:length(geol)
        if str2double(rmseg{i,1})==geol(k).FaultID
            lontemp = geol(k).lonseg;
            lattemp = geol(k).latseg;
            hp1=plot(lontemp,lattemp,'g-','linewidth',2);
%             hp2=plot(lontemp(rmseg{i,2}),lattemp(rmseg{i,2}),'r-','linewidth',2);

%             lontemp(rmseg{i,2})=[];
%             lattemp(rmseg{i,2})=[];
%             geol(k).lonseg=lontemp;
%             geol(k).latseg=lattemp;

%             pause;

        end
    end
end

for i=1:Rf
    for k=1:length(geol)
        if str2double(rmseg{i,1})==geol(k).FaultID
            lontemp = geol(k).lonseg;
            lattemp = geol(k).latseg;
            lontemp(rmseg{i,2})=[];
            lattemp(rmseg{i,2})=[];
            geol(k).lonseg=lontemp;
            geol(k).latseg=lattemp;
        end
    end
end

% Remove these entire faults
% they are mostly small and west of the SAF in socal
% Newport Inglewood, pita point, San Luis Bay, Yorba Linda (crosses Whittier) , San Mateo
% Santa Susana east,connector, Puenta Hills, Pitas point, Hector Mine, Big
% Pine , Compton, Carlsbad

rmflt = {'183','213','250','303','340','265','219','220','221','211','310','333','113','16','43','29','195'};
kout=[];
for i=1:length(rmflt)
    for k=1:length(geol)
        if str2double(rmflt{i})==geol(k).FaultID
            kout = [kout;k];
            lontemp = geol(k).lonseg;
            lattemp = geol(k).latseg;
            hp1=plot(lontemp,lattemp,'c-','linewidth',2);
        end
    end
end
geol(kout)=[];


%% find fault segments that intersect one another

disp(' ');
disp('Checking for fault segments that cross...');

crosslist=[];

cnt = 0;

for i=1:length(geol)
    lons = geol(i).lonseg;
    lats = geol(i).latseg;
    for j=1:(length(lons)-1)
        for q=1:length(geol)
            if q<i
                lonsc = geol(q).lonseg;
                latsc = geol(q).latseg;
                for z=1:(length(lonsc)-1)
                    [ido,lati,loni]=dotheyintersect(lats(j),lons(j),lats(j+1),lons(j+1),latsc(z),lonsc(z),latsc(z+1),lonsc(z+1));
                    if ido==1
                        %[Dist,~]=distance(lati,loni,[lats(j);lats(j+1);latsc(z);latsc(z+1)],[lons(j);lons(j+1);lonsc(z);lonsc(z+1)]);
                        %  if min(Dist)==0 && mamax(Dist)>0   % catches case where faults are end to end
                        %  if mindist ~= 0 &
                        disp(['Warning: geol(' num2str(i) ') segment ' num2str(j) ' intersects geol(' num2str(q) ') segment ' num2str(z)]);
                        crosslist=[crosslist; i j q z];
                        
                        figure(20);
                        clf;
                        plot(lons,lats,'g-o','linewidth',2);
                        text(median(lons),median(lats),num2str(i),'color','g');
                        text(lons(1),lats(1),'beg','color','g')
                        text(lons(end),lats(end),'end','color','g')
                        hold on;
                        plot(lonsc,latsc,'r-o','linewidth',2);
                        plot(loni,lati,'kp','markersize',12);
                        text(median(lonsc),median(latsc),num2str(q),'color','r');
                        text(lonsc(1),latsc(1),'beg','color','r')
                        text(lonsc(end),latsc(end),'end','color','r')
%                         pause;
                        cnt=cnt+1;
                    end
                end
            end
        end
    end
end

if cnt==0
    disp(' ');
    disp('No crossing faults!');
end

    
%%

save('./AfterBuildSources','geol');



