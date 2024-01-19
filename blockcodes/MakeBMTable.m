function MakeBMTable(file,slip,s_slip,omega,s_omega,creep,blocks,faults,nodes,lon,lat)
% MakeBMTable(file,slip,s_slip,omega,s_omega,creep,blocks,faults,nodes,lon,lat)
%
% Makes a text file with two tables.  One for block rotation rates.
% Another for fault slip rates. Two columns are given for fault-perpendicular slip rates,
% one for dip slip rates and anther for horizontal extension (where dip
% slip rate has been projected to the horizontal according to the fault dip).


eval(['!rm -rf ' file]);

fid=fopen(file,'w');

[ns,~]=size(slip);

if isempty(creep)
   fprintf(fid,'%s\n','Fault  StrikeSlip      DipSlip          Extension'); 
else
   fprintf(fid,'%s\n','Fault  StrikeSlip      Creep  Total   DipSlip        Extension'); 
end

for i=1:ns
    
    if ~isempty(creep)
        
        j=find(i==creep(:,1));
        
        if isempty(j)
            crpval = 0;
        else
            crpval = creep(j,2)*1000;
        end
        
        str = [sprintf('%4s',num2str(i)) '   ' ...
            sprintf('%6.2f',slip(i,1)*1000) '+/-' sprintf('%4.2f',s_slip(i,1)*1000) ...
            '   ' ...
            sprintf('%6.2f',crpval)  ...
            '   ' ...
            sprintf('%6.2f',slip(i,1)*1000 + crpval)  ...
            '   ' ...
            sprintf('%6.2f',slip(i,2)*1000) '+/-' sprintf('%4.2f',s_slip(i,2)*1000) ...
            '    ' ...
            sprintf('%6.2f',slip(i,2)*1000*cosd(faults(i,4))) '+/-' sprintf('%4.2f',s_slip(i,2)*1000*abs(cosd(faults(i,4))))];
        
    else
        str = [sprintf('%4s',num2str(i)) '   ' ...
            sprintf('%6.2f',slip(i,1)*1000) '+/-' sprintf('%4.2f',s_slip(i,1)*1000) ...
            '   ' ...
            sprintf('%6.2f',slip(i,2)*1000) '+/-' sprintf('%4.2f',s_slip(i,2)*1000) ...
            '    ' ...
            sprintf('%6.2f',slip(i,2)*1000*cosd(faults(i,4))) '+/-' sprintf('%4.2f',s_slip(i,2)*1000*abs(cosd(faults(i,4))))];
        
    end
    
    fprintf(fid,'%s\n',str);
    
end

fprintf(fid,'%s\n',' ');

[no,~]=size(omega);

[N]=SitesPerBlock(blocks,nodes,lon,lat);

R = 6371000;

fns = fieldnames(blocks);

fprintf(fid,'%s\n','Block  PoleLat         PoleLon      WRate(deg/Myr)   LatLonCorr    N  VertAxis(deg/Myr)'); 
for i=1:no
    
    wrate = norm(omega(i,:));
    swrate = sqrt((s_omega(i,1).^2)*(omega(i,1)/wrate).^2 + ...
        (s_omega(i,2).^2)*(omega(i,2)/wrate).^2 + ...
        (s_omega(i,3).^2)*(omega(i,3)/wrate).^2);
    
    OM = R*omega(i,:)/wrate;
    SOM = R*s_omega(i,:)/wrate;
    [lat,lon,~,slon,slat,~,covL]=xyz2latlon(OM(1),OM(2),OM(3),SOM(1),SOM(2),SOM(3),0,0,0,0);
    corr = covL(1,2)*(180/pi).^2/(slon*slat);

    %N=length(getfield(blocks,char(fns(i))));

    
    covL =  diag([1 1 1]);
    eval(['inodes=blocks.' char(fns(i)) ';']);
    
    latmid = mean(nodes(inodes,2));
    lonmid = mean(nodes(inodes,1));
    
    [x,y,z,~]=latlon2xyz(latmid,lonmid,0,covL);
    r = [x;y;z];
    r = r./norm(r);
    omega_vertrate = dot(omega(i,:),r);
    
    str = [sprintf('%4s',num2str(i)) '  ' sprintf('%5.1f',lat) '+/-' sprintf('%-6.1f',slat) ...
        '  ' sprintf('%6.1f',lon) '+/-' sprintf('%-6.1f',slon) '  ' ...
        sprintf('%5.3f',wrate*1e6*180/pi) '+/-' sprintf('%5.3f',swrate*1e6*180/pi) ...
        '  ' sprintf('%6.3f',corr) '     ' sprintf('%3.0f',N(i)) '   ' sprintf('%6.3f',omega_vertrate*1e6*180/pi)];
   
    fprintf(fid,'%s\n',str);
    
end



fclose(fid);




