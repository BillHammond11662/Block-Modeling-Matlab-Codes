function [nodes,blocks,faults,geol,ifltdamp]=MakeBlockModel(geol,bnds,params,optverb,optplot,optdamp)
% [nodes,blocks,faults,geol,ifltdamp]=MakeBlockModel(geol,bnds,params,optverb,optplot,optdamp)
%
% geol is the databae of fault properties and geometry
% bnds is the 4x1 boundary box
% optverb = 1 gives informative output
% optplot = 1 makes many plots of steps
% optdamp = 1 is used to add fault paramters for every block boundary (not
%            just the boundaries that represent real geologic structures in the fault database
% params is a cell array containing:
%   Bnmax        - maximum number of blocks 
%   IntAngleMax  - maximum allowed interior block polygon angle
%   IntAngleMin  - minimum allowed interior block polygon angle
%   Amin = 400   - minimum block area in km^2
%   Rmax = 1.9   - maximum shape parameter allowed
%   Degsmall=1   - minimum interior angle for sliver blocks on boundary
%


ifltdamp=[];

% trunctate the list of segments to those touching the domain.

lonp=[];
latp=[];
ikp = [];
for i=1:length(geol)
    
    lons = geol(i).lonseg;
    lats = geol(i).latseg;
    B=length(lons);
    keep=0;
    for j=1:(B-1)
        lon1=lons(j);
        lat1=lats(j);
        lon2=lons(j+1);
        lat2=lats(j+1);
        
        if ((lon1>bnds(1) && lon1<bnds(2) && lat1>bnds(3) && lat1<bnds(4)) || ...
            (lon2>bnds(1) && lon2<bnds(2) && lat2>bnds(3) && lat2<bnds(4)))
            keep=1;
        end
    end
    if keep==1
        ikp = [ikp;i];
        
        lonad = geol(i).lonseg;
        latad = geol(i).latseg;
        if size(lonad,2)>1 && size(lonad,1)==1
           lonad=lonad';
        end
        if size(latad,2)>1 && size(latad,1)==1
           latad=latad';
        end
        if size(lonad,2)>1 && size(lonad,1)>1
            error('Problem with geol.lonseg');
        end
        lonp = [lonp;lonad;nan];
        latp = [latp;latad;nan];
    end
end

geol=geol(ikp);

%% remove duplicate nodes in lonseg,latseg

for i=1:length(geol)
    lons = geol(i).lonseg;
    lats = geol(i).latseg;
    
    jout=[];
    for j=2:length(lons)
       if (lons(j)==lons(j-1) && lats(j)==lats(j-1))
           jout = [jout;j];
       end
    end
    lons(jout)=[];
    lats(jout)=[];
    geol(i).lonseg = lons;
    geol(i).latseg = lats;
end


%%  Clip the faults crossing outside the domain

% the faults sometimes cross out of our box so we'll need to trim the
% faults to fit in the box.

geol=ClipFaults(geol,bnds);


%% delete any fault trains that only have one point

iout=[];
for i=1:length(geol)
    if length(geol(i).lonseg)<=1
        iout=[iout;i];
    end
end

geol(iout)=[];

%% Find, warn about, repair fault segments that intersect one another

if optplot==1
    A=gcf;figure(A.Number+1);
    clf;
    for i=1:length(geol)
        
        lons = geol(i).lonseg;
        lats = geol(i).latseg;
        
        plot(lons,lats,'ro-','linewidth',2)
        hold on;
        text(lons(end),lats(end),num2str(i))
    end
    axis(bnds+[-.2 .2 -.2 .2]);
    set(gca,'dataaspectratio',[1 cosd(38) 1])
    xlabel('Longitude');
    ylabel('Latitude');
end


iter = 1;
still=1;
disp(' ');
while still
    disp(['Iterating to repair fault segments ' num2str(iter)]);
    
    crosslist=[];
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
                            [Dist,~]=distance(lati,loni,[lats(j);lats(j+1);latsc(z);latsc(z+1)],[lons(j);lons(j+1);lonsc(z);lonsc(z+1)]);
                            mindist = min(Dist);  % catches case where faults are end to end
                            if mindist ~= 0
                                if optverb==1
                                    disp(['Warning: geol(' num2str(i) ') intersects geol(' num2str(q) ')']);
                                    disp('Unpredictable or illegal blocks could form.');
                                    disp('Deciding which segments to remove.');
                                end
                                crosslist=[crosslist; i j q z];
                            end
                        end
                    end
                end
            end
        end
    end
    
    if isempty(crosslist)
        still=0;
    else
        
        % find which segment appears in the list most
        clu = unique(crosslist(:,1:2),'rows');
        jcnt=zeros(size(clu,1),1);
        for j=1:size(clu,1)
            jcnt(j) = length(find(clu(j,1)==crosslist(:,1) & clu(j,2)==crosslist(:,2)));
        end
        imax = find(jcnt==max(jcnt));
        fout=crosslist(imax(1),1);
        sout=crosslist(imax(1),2);
        
        if optverb
           disp(['Removing segment ' num2str(sout) ' from fault ' num2str(fout)]);
        end
        % remove this segment
        %    if segment is end or beginning, delete the one node needed to shorten
        %    if segment is in middle truncate geol(i) to geol(1:j) and add new geol(j+1:end)
        
        lons = geol(fout).lonseg;
        lats = geol(fout).latseg;
        SL=length(lons);
        if sout==1
            lons(1)=[];
            lats(1)=[];
            geol(fout).lonseg=lons;
            geol(fout).latseg=lats;
        elseif sout==(SL-1)
            lons(end)=[];
            lats(end)=[];
            geol(fout).lonseg=lons;
            geol(fout).latseg=lats;
        else
            lons2=lons(sout+1:end);
            lats2=lats(sout+1:end);
            lons(sout+1:end)=[];
            lats(sout+1:end)=[];
            geol(fout).lonseg=lons;
            geol(fout).latseg=lats;
            geol(length(geol)+1)=geol(fout);
            geol(length(geol)).lonseg=lons2;
            geol(length(geol)).latseg=lats2;
        end
        
        if optplot==1
            A=gcf;figure(A.Number+1);
            clf;
            for i=1:length(geol)
                
                lons = geol(i).lonseg;
                lats = geol(i).latseg;
                
                plot(lons,lats,'ro-','linewidth',2)
                hold on;
                text(lons(end),lats(end),num2str(i))
            end
        end
        
    end
    iter = iter+1;
end



%% create evenly spaced nodes along faults to structure block boundaries 

%dint = .1; % in degrees
% dint = .05; % in degrees
dint = min([bnds(2)-bnds(1) bnds(4)-bnds(3)])/40;

for i=1:length(geol)
    
    [D,~]=distance( geol(i).latseg(1), geol(i).lonseg(1),geol(i).latseg(end), geol(i).lonseg(end));
    c=ceil(D/dint);
    
    m=length(geol(i).latseg);
    
    latint=interp1(1:m,geol(i).latseg,linspace(1,m,c));
    lonint=interp1(1:m,geol(i).lonseg,linspace(1,m,c));
    
    if length(latint)>1
        geol(i).latseg=latint';
        geol(i).lonseg=lonint';
    end
end

if optplot==1
    A=gcf;figure(A.Number+1);
    clf;
    for i=1:length(geol)
        
        lons = geol(i).lonseg;
        lats = geol(i).latseg;
        
        plot(lons,lats,'ro-','linewidth',2)
        hold on;
        text(lons(end),lats(end),num2str(i))
    end
end


%% Plot the resulting fault segments

lonps=[];
latps=[];
lonpsplot=[];
latpsplot=[];
for i=1:length(geol)
   lonps = [lonps;geol(i).lonseg];
   latps = [latps;geol(i).latseg];
   lonpsplot = [lonpsplot;NaN;geol(i).lonseg];
   latpsplot = [latpsplot;NaN;geol(i).latseg];
end

lonbnd = [bnds(1) bnds(2) bnds(2) bnds(1) bnds(1)];
latbnd = [bnds(3) bnds(3) bnds(4) bnds(4) bnds(3)];

if optplot==1
    A=gcf;figure(A.Number+1);
    clf;
    
    PlotStates(0,0);
    hold on;
    plot(lonp,latp,'-','linewidth',2,'color',[.2 .2 .5]);
    %plotm(latps,lonps,'ko','markersize',12,'markerfacecolor','r');
    plot(lonpsplot,latpsplot,'ko','markersize',10,'markerfacecolor','r');
    plot(lonbnd,latbnd,'k-');
    axis(bnds+[-.2 .2 -.2 .2]);
    set(gca,'dataaspectratio',[1 cosd(38) 1])
    xlabel('Longitude');
    ylabel('Latitude');
end


%% add points around perimeter

n=4; % number of segments along which to divide each boundary
lonadd = linspace(bnds(1),bnds(2),n+1)'; 
for i=2:n
    lonadd=[lonadd;bnds(1);bnds(2)];
end
lonadd = [lonadd;linspace(bnds(1),bnds(2),n+1)'];

latvals = linspace(bnds(3),bnds(4),n+1);
latadd(1:n+1,1)=bnds(3);
for i=2:n
    latadd=[latadd;latvals(i);latvals(i)];
end
latadd = [latadd;bnds(4)*ones(n+1,1)];

% remove elements of lonadd that are too close to lonps

% keyboard;

% l1=length(latadd);
% l2=length(lonps);
% D=nan(l1,l2);
% dmin = min([bnds(2)-bnds(1);bnds(4)-bnds(3)])/50;  % in degrees
% iout=[];
% for i=1:l1
%     [D,~]=distance(latadd(i),lonadd(i),latps,lonps);
%     if (min(D)<dmin)
%         iout=[iout;i];
%     end
% end
% latadd(iout)=[];
% lonadd(iout)=[];

lonps=[lonps;lonadd];
latps=[latps;latadd];


%% Delaunay triangulation

nodes=[lonps latps];
nodes=unique(nodes,'rows');
lonps=nodes(:,1);
latps=nodes(:,2);

tri = delaunay(lonps,latps);
[T,~]=size(tri);

%% create triangulation segment list from triangle list and plot it

lonlatseg=[];
for i=1:T
    
    n1=tri(i,1);
    n2=tri(i,2);
    n3=tri(i,3);
    
    for j=1:3
        nn = [j;j+1];
        nn = mod(nn,3);
        nn(nn==0)=3;
        np = tri(i,nn);
        
        %only add the segment if it is not already present in the list 
        if i==1
            isit = [];
        else
            isit = find((lonlatseg(:,1)==lonps(np(1)) & lonlatseg(:,2)==latps(np(1))  &  lonlatseg(:,3)==lonps(np(2)) & lonlatseg(:,4)==latps(np(2)) ) | ...
                        (lonlatseg(:,3)==lonps(np(1)) & lonlatseg(:,4)==latps(np(1))  &  lonlatseg(:,1)==lonps(np(2)) & lonlatseg(:,2)==latps(np(2)) ));
        end
        
        if isempty(isit)
            lonlatseg=[lonlatseg;lonps(np(1)) latps(np(1)) lonps(np(2)) latps(np(2))];
        end
    end
end

if optplot==1
    A=gcf;figure(A.Number+1);
    clf;
    
    plot(lonp,latp,'-','linewidth',2,'color',[.2 .2 .5]);
    hold on;

    %plotm(latps,lonps,'ko','markersize',12,'markerfacecolor','r');
    plot(lonpsplot,latpsplot,'ko','markersize',10,'markerfacecolor','r');
    plot(lonbnd,latbnd,'k-');

    for i=1:length(lonlatseg)
        plot(lonlatseg(i,[1 3]),lonlatseg(i,[2 4]),'m-');
         text(mean(lonlatseg(i,[1 3])), mean(lonlatseg(i,[2 4])),num2str(i));
    end
    set(gca,'dataaspectratio',[1 cosd(38) 1])
    xlabel('Longitude');
    ylabel('Latitude');
    axis(bnds);
end

% keyboard;


%% if any fault segment is crossed by a triangle side, then
% divide that fault segment in two and redo the triangulation.
% 
% 
% if optplot==1
%     A=gcf;figure(A.Number+1);
%     clf;
%     hold on;
%     plot(lonp,latp,'-','linewidth',2,'color',[0 0 1]);
%     plot(lonps,latps,'ko','markersize',10,'markerfacecolor','r');
%     axis(bnds);
%     
%     trimesh(tri,lonps,latps,'color',[.5 .5 .5]);
%     
%     set(gca,'dataaspectratio',[1 cosd(38) 1])
%     xlabel('Longitude');
%     ylabel('Latitude');
% end


%% make the first block model
% it has a lot of small triangular blocks

[G,~]=size(lonlatseg);

% first make a block model utilizing ALL triangles in Delaunay Tri.
[M,~]=size(tri);
blocks=[];
for i=1:M
  blocks.(['b' num2str(i)]) = [tri(i,1) tri(i,2) tri(i,3) tri(i,1)];
end

% next make fault list, adding only those that conform to fault segments
faults =[];
for i=1:G
    lon1=lonlatseg(i,1);lat1=lonlatseg(i,2);
    lon2=lonlatseg(i,3);lat2=lonlatseg(i,4);
    
    iflt=seg2flt(lon1,lat1,lon2,lat2,geol);
    ison=isonboundary(lat1,lon1,lat2,lon2,bnds);

    if ~isempty(iflt) && ~ison
        
        n1 = find(lon1==nodes(:,1) & lat1==nodes(:,2));
        n2 = find(lon2==nodes(:,1) & lat2==nodes(:,2));
        if isempty(n1) || isempty(n2) || ischar(n1) || ischar(n2) 
            keyboard;
        end
        
        % remember:
        %  positive fault dip is down to the left while walking from n1 to n2.
        %  negative fault dip is down to the right while walking from n1 to n2.
        
        dip = geol(iflt).dip;
        if ischar(dip)
            dip = str2double(dip);
        end
        if dip==0
            diptmp = char(geol(iflt).disp_dips);
            dip = str2double(diptmp(1:2));
        end
        dipdir = uppercase(geol(iflt).dip_dir);
        if strcmp('E',dipdir)
            dipdirdd = 90;
        elseif strcmp('W',dipdir)
            dipdirdd = 270;
        elseif strcmp('N',dipdir)
            dipdirdd = 0;
        elseif strcmp('S',dipdir)
            dipdirdd = 180;
        elseif strcmp('NE',dipdir)
            dipdirdd = 45;
        elseif strcmp('SE',dipdir)
            dipdirdd = 135;
        elseif strcmp('NW',dipdir)
            dipdirdd = 315;
        elseif strcmp('SW',dipdir)
            dipdirdd = 225;
        elseif strcmp('Vertical',dipdir) || strcmp('VERTICAL',dipdir)
            dipdirdd = [];  % this will force dip to 90˚
        else
            dipdirdd = [];
        end
        
        [~,az]=distance(nodes(n1,2),nodes(n1,1),nodes(n2,2),nodes(n2,1));
        
        if ~isempty(dipdirdd)
            angdif = mod(az-dipdirdd,360);
            if (180<angdif && angdif<=360)  % this happens when the dip is to the right
                dip=-1*dip;
            end
        else
            dip=90;
        end
        
        if isempty(dip)
            keyboard;
        end
        
        faults=[faults;
            n1 n2 15 dip];

    end
end

if optplot==1
    
    A=gcf;figure(A.Number+1);
    clf
    PlotBlocks(nodes,blocks,faults,1,1,0);
    axis(bnds);
    title([num2str(length(fieldnames(blocks))) ' Blocks']);
    
    A=gcf;figure(A.Number+1);
    clf
    PlotFaults(blocks,nodes,faults,0,0,0);
    axis(bnds);
    title([num2str(length(fieldnames(blocks))) ' Blocks, ' num2str(size(faults,1)) ' Faults']);
end

%%  perform imporant checks on model

[~,four] = size(faults);
if four~=4
    %keyboard;
    blocks=[];
    ifltdamp=[];
    disp('Warning. Null model.  No blocks or faults.');
    return;
end

if ~isempty(faults)
    flag=CheckFaults(nodes,blocks,faults);
    if flag==1
        disp('Block model developed illegal faults.. Aborting construction...');
        blocks=[];
        ifltdamp=[];
        return;
    end
end

CheckBlocks(blocks);

%%

[score,lonlatseg]=ScoreModel(blocks,nodes,faults,bnds);

if optplot==1
    A=gcf;figure(A.Number+1);
    clf
    PlotScore;
end


%% Iterate

% Algorithm to find block model with
% 1) not too many blocks
% 2) few very small blocks
% 3) few blocks with small aspect ratios (shape parameters)
% 4) no blocks with very large internal angles
% 5) no blocks with very small internal angles

Bnmax = params{1};
IntAngleMax = params{2};
IntAngleMin = params{3};
Amin = params{4};
Rmax = params{5};
Degsmall = params{6};

still = 1;
iter = 1;
while still
    
    disp(' ');
    disp(['Iterating to combine blocks to simplify model: ' num2str(iter)]);
    disp(' ');

    % remove segements with score 0, as long as its removal does not
    % generate a block with a large internal angle
    if (optverb==1)
        disp(' ');
        disp('Removing blocks by score.');
    end
    
    if iter>1
        [score,lonlatseg]=ScoreModel(blocks,nodes,faults,bnds);
    end
        
    % the first round removes block boundaries by score
    %
    if any(find(score==0))
        
        [blocks,faults]=ReduceBlockByScore(blocks,faults,nodes,lonlatseg,score,optverb);
        
        if isempty(blocks)
            return;
        end
        
        if optplot==1
            
            [score,lonlatseg]=ScoreModel(blocks,nodes,faults,bnds);
            A=gcf;figure(A.Number+1);
            PlotScore;
            A=gcf;figure(A.Number+1);
            PlotBlockIter;
        end
    else
        disp(' ');
        disp('No need to remove any with score==0.')
    end
        
    % identify and remove sliver blocks on boundary
    % if they are on a boundary, just delete them (don't combine them).
    if (optverb==1)
        disp(' ');
        disp('Removing sliver boundary blocks.');
    end
    [blocks,faults]=RemoveSliverBlocks(blocks,faults,nodes,bnds,Degsmall,optverb);
    if optplot==1
        A=gcf;figure(A.Number+1);
        PlotBlockIter;
    end
    
    % identify and remove segments associated with small blocks, as long as its removal does not
    % generate a block with a large internal angle
    if (optverb==1)
        disp(' ');
        disp('Removing blocks by size.');
    end
    [blocks,nodes,faults]=ReduceBlockBySize(blocks,nodes,faults,Amin,bnds,optverb);
    
    if optplot==1
        A=gcf;figure(A.Number+1);
        PlotBlockIter;
    end
            
    % identify and remove segments associated with interior angles outside allowed range
    if (optverb==1)
        disp(' ');
        disp('Removing blocks by interior angles.');
    end
    [blocks,faults]=ReduceBlockByInteriorAngles(blocks,nodes,faults,bnds,IntAngleMin,IntAngleMax,optverb);
    if optplot==1
        A=gcf;figure(A.Number+1);
        PlotBlockIter;
    end
    
    % identify and remove segments associated with and small aspect ratio blocks, as long as its removal does not
    % generate a block with a large internal angle
    if (optverb==1)
        disp(' ');
        disp('Removing blocks by Shape Parameter.');
    end
    [blocks,nodes,faults]=ReduceBlockByShapeParameter(blocks,nodes,faults,bnds,Rmax,IntAngleMin,IntAngleMax,optverb);
    if optplot==1
        A=gcf;figure(A.Number+1);
        PlotBlockIter;
    end
    
    Mlast=M;
    M=length(fieldnames(blocks));
    if (M<Bnmax)
        still=0;
        disp(' ');
        disp('Number of Blocks small enough. Quitting this nonsense!');
        disp(' ');
    end
    if Mlast==M
        still=0;
        disp(' ');
        disp('Number of Blocks not changing. Quittin Time!');
        disp(' ');
    end
    if size(faults,1)==0
        still=0;
        disp(' ');
        disp('Uh oh. All faults removed. No point now. Stopping');
        disp(' ');
        blocks=[];
        return;
    end
    
    iter=iter+1;
    
    % Block names can get a little long when they are defined by the
    % combination of two blocks.  So here we rename the blocks with neater short
    % names.
    blocks = NiceBlockNames(blocks);
    
end

[~,~]=CheckNodes(nodes,faults,blocks);


%%  delete any fault segments that lie on the boundary, or have less than two adjacent blocks 

if optverb
    disp(' ');
    disp('Deleting Faults that are on the boundary');
end

P=size(faults,1);
fout=[];
for i=1:P
    j=Fault2Blocks(nodes,blocks,faults,i);
    if length(j)<2
        fout =[fout;i];
        disp(['Deleting block ' num2str(i)]);
    end
end
if ~isempty(fout)
    faults(fout,:)=[];
end

%% Look for faults that are 'crossed' by another through-going boundary that is not necessarily a fault
%  this prevents splitting faults like the Garlock with a NW trending block
%  boundary.
%  do this by:
%     scanning each fault
%     looking for nodes that are associated with 4 or more blocks
%     and combining blocks that don't result in eliminating faults

blocknames=fieldnames(blocks);
M=length(blocknames);
segrmlist=[];
ndlist=[];
[~,ndseg]=Model2segs(blocks,nodes);
[A,~,~,~,~,~,~,~] = blockgeoms(blocks,nodes);

for i=1:length(geol)

    lonf = geol(i).lonseg;
    latf = geol(i).latseg;

    if length(lonf)>2
        for j=2:(length(lonf)-1)

            nd = find(nodes(:,1)==lonf(j) & nodes(:,2)==latf(j),1);

            % make sure its part of a fault and not on end
            

            if ~isempty(nd)
                bnum=[];
                for k=1:M
%                     bnodes = eval(['blocks.' blocknames{k}]);
                    bnodes = blocks.(blocknames{k});
                    if ~isempty(intersect(nd,bnodes))
                        bnum=[bnum;k];
                    end
                end
                bnum=unique(bnum);
                if length(bnum)>=4

                    ndlist = [ndlist;nd];

%                     keyboard;

                    % find all segments connected to this node nd
                    [nseg,~]=find(ndseg==nd);

                    % go through each segment connected to the node see if removing it is an
                    % option. Only consider segments that are not faults in
                    % the geol database for removal
                    segrmlist0=[];
                    for z=1:size(nseg,1)

                        n1=ndseg(nseg(z),1);
                        n2=ndseg(nseg(z),2);
                        iflt=seg2flt(nodes(n1,1),nodes(n1,2),nodes(n2,1),nodes(n2,2),geol);

                        if isempty(iflt)
                            % now add to the list of segments to
                            % potentially remove
                            segrmlist0=[segrmlist0;nseg(z)];
                        end
                    end

%                     keyboard;

                    % for this node go through remaining possible segments and find best to remove.
                    % It is the one with a legal combination (fout not empty), and with the minumum area
                    % of the combined blocks 
                    Acomb = nan(length(segrmlist0),1);
                    for y = 1:length(segrmlist0)
                        segnds = ndseg(segrmlist0(y),:)';  % these are the nodes of the segment
                        [~,~,bout,~]=RemoveBlock(blocks,nodes,faults,segnds);
                        if ~isempty(bout)
                           Acomb(y,1)=sum(A(bout));
                        end
                    end
                    [~,iminA]=min(Acomb);  % if Acomb is not all NaN then iminA is the least Area combined block
                    if ~isnan(iminA)  
                        segrmlist=[segrmlist;segrmlist0(iminA)];
                    end
                    
                end
            end
        end
    end
end

% hold on;
% plot(nodes(ndlist,1),nodes(ndlist,2),'bo','markerfacecolor','b','markersize',12);
% for i=1:length(segrmlist)
%     n1=ndseg(segrmlist(i),1);
%     n2=ndseg(segrmlist(i),2);
%     plot(nodes([n1 n2],1),nodes([n1 n2],2),'b-','linewidth',2)
% end

if optverb
    disp(' ');
    disp(['Removing ' num2str(length(segrmlist)) ' segments to aviod block boundaries crossing geol faults.']);
end

% keyboard;

% Now actually remove those segments by combining blocks
for i=1:length(segrmlist)
    segnds=ndseg(segrmlist(i),:)';
    [blocksnew,faultsnew,bout,~]=RemoveBlock(blocks,nodes,faults,segnds);
    if ~isempty(bout)
        blocks = blocksnew;
        faults = faultsnew;
    end
end

if optplot==1
    A=gcf;figure(A.Number+1);
    PlotBlockIter;
end

    

%%  add faults to interior block boundaries and remember which ones they are
%  create an option for this.

if optdamp == 1
    Nf = size(faults,1);
    bnames = fieldnames(blocks);
    M=length(bnames);
    
    for i=1:M
        
        %iterate through each boundary
%         eval(['jj=blocks.' char(bnames(i)) ';']);
        jj=blocks.(bnames{i});
        for k = 2:length(jj)
            n1=jj(k);
            n2=jj(k-1);
            
            % is added to the fault list yet?
            iflt=nodes2fault(faults,n1,n2);
            % are there two blocks?
            blist=seg2blocks(blocks,[n1;n2]);

            % length of fault in km
            flen = deg2km(baz(nodes(n1,2),nodes(n1,1),nodes(n2,2),nodes(n2,1)));
            % segment needs to be longer than 100 meters.
            if isempty(iflt) && length(blist)>=2 && flen>0.100 
                faults = [faults;
                          n1 n2 15 45];
            end
        end
    end
    ifltdamp = (Nf+1):size(faults,1);
end



%% Relabel blocks with easier names

blocks = NiceBlockNames(blocks);

% bnamesold = fieldnames(blocks);
% for i=1:length(fieldnames(blocks))
%    bnamenew = ['block' sprintf('%.0f',i)];
%    ndlist = blocks.(bnamesold{i});
%    blocks.(bnamenew)=ndlist;
%    blocks=rmfield(blocks,bnamesold{i});
% end

if optplot==1
    A=gcf;figure(A.Number+1);
    PlotBlockIter;
end

[flag,isus]=CheckNodes(nodes,faults,blocks);


%% Check faults again to avoid NaNs.
%
% sometime later need to go back and find where Nans in faults are
% originating from

if ~isempty(faults)
    flag=CheckFaults(nodes,blocks,faults);
    if flag==1
        disp(' ');
        disp('Block model developed illegal faults.. Aborting construction...');
        blocks=[];
        ifltdamp=[];
        faults=[];
        return;
    end
end






