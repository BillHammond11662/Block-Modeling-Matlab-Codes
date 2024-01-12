function omegavaxis=SpinRate(blocks,nodes,omega)
% omegavaxis=SpinRate(blocks,nodes,omega)
%
% rotation rate for each block in radians/yr


blocknames =fieldnames(blocks);

for i=1:length(blocknames)

    ndlist=getfield(blocks,char(blocknames(i)));
    lon0=mean(nodes(ndlist,1));
    lat0=mean(nodes(ndlist,2));
    [rx,ry,rz,dum]=latlon2xyz(lat0,lon0,0,[]);
    R=norm([rx ry rz]);
    rx=rx./R;
    ry=ry./R;
    rz=rz./R;
    
    omegavaxis(i,1)=  omega(i,1)*rx + omega(i,2)*ry + omega(i,3)*rz;

end