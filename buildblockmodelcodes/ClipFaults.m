function geol=ClipFaults(geol,bnds)
% geol=ClipFaults(geol,bnds)
%
% 

lonbnd = [bnds(1) bnds(2) bnds(2) bnds(1) bnds(1)];
latbnd = [bnds(3) bnds(3) bnds(4) bnds(4) bnds(3)];

for i=1:length(geol)
    
    lons = geol(i).lonseg;
    lats = geol(i).latseg;
    decision=zeros(size(lons));
    
    % 0 delete
    % 1 keep
    % 2 move (but also keep)
    decision(inpolygon(lons,lats,lonbnd,latbnd))=1;
    
    % figure out which segment nodes are outside but connected to node inside
    for j=1:length(lons)
        if decision(j)~=1
            if j==1
                if inpolygon(lons(j+1),lats(j+1),lonbnd,latbnd) 
                   decision(j)=2;
                end
            elseif j<length(lons)
                if (inpolygon(lons(j+1),lats(j+1),lonbnd,latbnd) || inpolygon(lons(j-1),lats(j-1),lonbnd,latbnd))
                    decision(j)=2;
                end
            else % j==length(lons)
                if inpolygon(lons(j-1),lats(j-1),lonbnd,latbnd)
                    decision(j)=2;
                    % move it
                end
            end
        end
    end
    
    lons(decision==0)=[];
    lats(decision==0)=[];
    decision(decision==0)=[];
    
    % move the nodes that are outside but connected to inside nodes
    lontmp=[];
    lattmp=[];
    split = 0;
    for j=1:length(lons)
        
        % if its in add this node to the list
        if decision(j)==1
            lontmp=[lontmp;lons(j)];
            lattmp=[lattmp;lats(j)];
            
        elseif decision(j)==2
            found = 0;
            
            % if its out, but next is in, move this one to boundary.
            if j<length(lons)
                if inpolygon(lons(j+1),lats(j+1),lonbnd,latbnd)
                    found=1;
                    [latint,lonint]=intersectbounds(lats(j),lons(j),lats(j+1),lons(j+1),bnds);
                    if length(latint)~=1 || length(lonint)~=1
                        keyboard;
                    end
                    lontmp=[lontmp;lonint];
                    lattmp=[lattmp;latint];
                end
            end
            
            % if its out, and previous is in
            if j>1
                if inpolygon(lons(j-1),lats(j-1),lonbnd,latbnd)
                    found=1;
                    [latint,lonint]=intersectbounds(lats(j-1),lons(j-1),lats(j),lons(j),bnds);
                    lontmp=[lontmp;lonint];
                    lattmp=[lattmp;latint];
                end
            end
                        
            % if out, and connected to previous and next.  This sometimes
            % happens. Split the noe
            if j>1 && j<length(lons)
                if (inpolygon(lons(j+1),lats(j+1),lonbnd,latbnd) && inpolygon(lons(j-1),lats(j-1),lonbnd,latbnd))
                    found=1;
                    [latint1,lonint1]=intersectbounds(lats(j-1),lons(j-1),lats(j),lons(j),bnds);
                    [latint2,lonint2]=intersectbounds(lats(j+1),lons(j+1),lats(j),lons(j),bnds);
                    lontmp=[lontmp;lonint1;lonint2];
                    lattmp=[lattmp;latint1;latint2];
                    split=1;
                end
            end
            
            if found==0
                error('decision==2 but not connected');
            end
        end
    end
    
    % remove duplicate nodes
    jout=[];
    for j=2:length(lontmp)
       if (lontmp(j)==lontmp(j-1) && lattmp(j)==lattmp(j-1))
           jout = [jout;j];
       end
    end
    lontmp(jout)=[];
    lattmp(jout)=[];

    geol(i).lonseg=lontmp;
    geol(i).latseg=lattmp;
    
%     if split==1
%         figure(90);
%         clf;
%         plot(lons,lats,'k-o');
%         hold on;
%         plot(bnds([1 2 2 1 1]),bnds([3 3 4 4 3]),'r-');
%         plot(lontmp,lattmp,'r*')
%         keyboard;
%     end
        
end


