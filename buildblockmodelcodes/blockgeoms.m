function [A,R,T,C,P,amax,amin,Nv] = blockgeoms(blocks,nodes)
% [A,R,T,C,P,amax,amin,Ns] = blockgeoms(blocks,nodes)
%   
% blocks, nodes are from block model
%
% A = area of each block in km^2
% R = sqrt(perimeter^2/area) - i.e. shape parameter. cirle=sqrt(4pi)~3.5
% T = Lmin/Lmax  - for triangles is aspect ratio
% C = 1 is convex, 0 is concave
% P = block perimeter
% amax,amin and max and min internal angles in the polygon
% Nv is number of vertices

% make histograms of block size and aspect ratio
bnames = fieldnames(blocks);
A=nan(size(bnames));
R=A;
T=A;
C=A;
amin=A;
amax=A;
Nv=A;
P=A;

nodesx=deg2km(nodes(:,1) - min(nodes(:,1)))*cosd(mean(nodes(:,2)));
nodesy=deg2km(nodes(:,2) - min(nodes(:,2)));

for i=1:length(bnames)
    nn=blocks.(bnames{i});

    if (length(nn)>2 && ~isempty(nn) && ~any(isnan(nodesx(nn))) && ~any(isnan(nodesy(nn))))
        
        blockshape  =polyshape(nodesx(nn),nodesy(nn),'KeepCollinearPoints',true);
        %[ GEOM, ~, ~ ] = polygeom( nodesx(nn), nodesy(nn));
        % GEOM = [ area   X_cen  Y_cen  perimeter ]
        % A(i) = GEOM(1);  
        A(i)=area(blockshape);
        % P(i) = GEOM(4);
        P(i) = perimeter(blockshape);
        % R(i) = sqrt(GEOM(4).^2/GEOM(1));
        R(i) = sqrt(P(i).^2/A(i));  % shape parameter 
        R(i) = R(i)/sqrt(4*pi);     % so a circle has R=1

        inn = unique(nn);
        if length(inn)==3  % then its a triangle
            n1=inn(1); n2=inn(2); n3=inn(3);
            lon1=nodes(n1,1);lat1=nodes(n1,2);
            lon2=nodes(n2,1);lat2=nodes(n2,2);
            lon3=nodes(n3,1);lat3=nodes(n3,2);
            
            [L1,~] = distance(lat1,lon1,lat2,lon2);
            [L2,~] = distance(lat2,lon2,lat3,lon3);
            [L3,~] = distance(lat3,lon3,lat1,lon1);
            
            Lmin = min([L1 L2 L3]);
            Lmax = max([L1 L2 L3]);
            T(i) = Lmin/Lmax;
            
        end
        
        C(i)= length(find(convexx([nodesx(nn) nodesy(nn)])==-1));
        if nn(1)==nn(end)
            if C(1)==C(end) && C(1)==-1
                C(i)=C(i)-1;
            end
        end
        
        angs = blockangles(nodes(nn,1),nodes(nn,2));
        amin(i) = min(angs);
        amax(i) = max(angs);
        
        Nv(i)=length(angs);
        
    end
end


