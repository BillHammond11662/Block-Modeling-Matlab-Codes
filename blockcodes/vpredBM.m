function [vepred,vnpred,vrote,vrotn,vstre,vstrn,vcose,vcosn]=...
          vpredBM(faults,nodes,blocks,L,nu,mu,omega,slip,strain,lon,lat)
%
% [vepred,vnpred,vrote,vrotn,vstre,vstrn,vcose,vcosn]=...
%         vpredBM(faults,nodes,blocks,L,nu,mu,omega,slip,strain,lon,lat);
%
% L is the number of fault segments allowed close to each point.
%
%

iln=find(lon>180);
lon(iln)=lon(iln)-360;

blocknames = fieldnames(blocks);
M = length(fieldnames(blocks));
Dtop =0;

R = 6371009.3;

% predicted velocities
plon = (nodes(faults(:,1),1) + nodes(faults(:,2),1))/2;
plat = (nodes(faults(:,1),2) + nodes(faults(:,2),2))/2;

vcosn=NaN(size(lon));
vcose=NaN(size(lon));
vrotn=NaN(size(lon));
vrote=NaN(size(lon));
vstrn=NaN(size(lon));
vstre=NaN(size(lon));

L = min([L size(faults,1)]);

for i=1:length(lon)

    bnum=FindBlock(blocks,nodes,lon(i),lat(i));

    if length(bnum)>1
        bnum=bnum(1);
    end

    if ~isempty(bnum)

        % Rotation Part
        [vrotn(i,1),vrote(i,1)]=PolePred(omega(bnum,:)',lat(i),lon(i));

        % Coseismic Part
        vcosn(i,1)=0;
        vcose(i,1)=0;
        dist=sqrt((lon(i)-plon).^2 + (lat(i)-plat).^2);
        [~,isrt]=sort(dist);

        if ~isempty(slip)
            for k=1:L

                k2 = isrt(k);
                Uss = [slip(k2,1)  0 0];
                Un  = [0  slip(k2,2) 0];
                lat1 = nodes(faults(k2,1),2);
                lon1 = nodes(faults(k2,1),1);
                lat2 = nodes(faults(k2,2),2);
                lon2 = nodes(faults(k2,2),1);
                dip  = faults(k2,4);
                W    = abs(faults(k2,3)/sind(dip));

                if ((lon(i)==lon1 && lat(i)==lat1) || (lon(i)==lon2 && lat(i)==lat2))
                    vlonn=0;vlonss=0;vlatn=0;vlatss=0;
                else
                    [vlonss,vlatss,vuss]=...
                        OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Uss,Dtop,dip,W,lon(i),lat(i),0);
                    [vlonn,vlatn,vun]=...
                        OkadaBlock(nu,mu,lat1,lon1,lat2,lon2,Un,Dtop,dip,W,lon(i),lat(i),0);
                    vlonn=real(vlonn);
                    vlonss=real(vlonss);
                    vlatn=real(vlatn);
                    vlatss=real(vlatss);
                end

                vcosn(i,1) = vcosn(i,1)+vlatss+vlatn;
                vcose(i,1) = vcose(i,1)+vlonss+vlonn;

            end
        end;

        % Strain Part

        if all(~isnan(strain(bnum,:)))

            % find centroid of block

            eval(['jj=blocks.' char(blocknames(bnum)) ';']);
            lat0 = mean(nodes(jj,2));
            lon0 = mean(nodes(jj,1));
            theta0 = deg2rad(90-lat0);

            dlon = deg2rad(lon(i) - lon0);  % in radians
            dlat = deg2rad(lat(i) - lat0);  % in radians
            dtheta=-dlat;

            e_phiphi = strain(bnum,1);
            e_thetaphi = strain(bnum,2);
            e_thetatheta = strain(bnum,3);

            % see Savage et al., 2001
            vstre(i,1) = R*dlon*e_phiphi*sin(theta0) + R*e_thetaphi*dtheta;
            vstrn(i,1) = -R*dlon*e_thetaphi*sin(theta0) - R*e_thetatheta*dtheta;
            
        else
            vstre(i,1) = 0;
            vstrn(i,1) = 0;
        end
    end
end

vnpred = vrotn-vcosn+vstrn;
vepred = vrote-vcose+vstre;
