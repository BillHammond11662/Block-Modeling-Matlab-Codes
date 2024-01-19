function [score,lonlatseg] = ScoreModel(blocks,nodes,faults,bnds)
% [score,lonlatseg] = ScoreModel(blocks,nodes,faults,bnds)
% 
%

lonlatseg=Model2segs(blocks,nodes);

%% make list of fault train end points

lonep=[];
latep=[];
P=size(faults,1);
for i=1:P
%     disp([num2str(i) ' of ' num2str(P)]);
    fnums=GetFaultsConnected(faults,i);
    nds = faults(fnums,1:2);
    ndsu = unique(nds);
    for k=1:length(ndsu)
        if length(find(ndsu(k)==nds))==1
            lonep=[lonep;nodes(ndsu(k),1)];
            latep=[latep;nodes(ndsu(k),2)];
        end
    end
end

[~,ia,~]=unique([lonep latep],'rows');
lonep=lonep(ia);
latep=latep(ia);

%%

G=size(lonlatseg,1);
score = zeros(G,1);

for i=1:G
    
    lon1=lonlatseg(i,1);
    lat1=lonlatseg(i,2);
    lon2=lonlatseg(i,3);
    lat2=lonlatseg(i,4);
    
    % if one end is touching a fault train endpoint
    if ~isempty(find(lon1==lonep & lat1==latep, 1))
      score(i)=score(i)+1;
    end 
    
    % if the other end is touching a fault train endpoint
    if ~isempty(find(lon2==lonep & lat2==latep, 1))
      score(i)=score(i)+1;
    end
        
    % if it is actually a fault segment
    if ~isempty(faults)
       iflt=seg2faultnum(lon1,lat1,lon2,lat2,faults,nodes);
    else
       iflt=[];
    end
    if ~isempty(iflt)
        score(i)=Inf;
    end
    
    % if its a boundary segment
    if ((lon1==bnds(1) && lon2==bnds(1)) || ...
        (lon1==bnds(2) && lon2==bnds(2)) || ...
        (lat1==bnds(3) && lat2==bnds(3)) || ...
        (lat1==bnds(4) && lat2==bnds(4)) )       
        score(i)=Inf;
    end

end

