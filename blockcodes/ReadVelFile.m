function [lat,lon,vn,ve,vu,sn,se,su,corrne,correu,corrnu,sta]=ReadVelFile(filename)
% function [lat,lon,vn,ve,vu,sn,se,su,corrne,correu,corrnu,sta]=ReadVelFile(filename)
%
% 

fid =fopen(filename);
C=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %s');
fclose(fid);

lat=C{1};
lon=C{2};
vn=C{3};
ve=C{4};
vu=C{5};
sn=C{6};
se=C{7};
su=C{8};
corrne=C{9};
correu=C{10};
corrnu=C{11};
sta=C{12};

