function Blocks2GMT(gmtfile,nodes,blocks)
% Blocks2GMT(gmtfile,nodes,blocks)

if exist(gmtfile,'file');
    delete(gmtfile);
end;

fid=fopen(gmtfile,'w');

blocknames = fieldnames(blocks);
[M,d]=size(blocknames);
for i=1:M

    eval(['jj=blocks.' char(blocknames(i)) ';']);

    for j=2:length(jj)

        lon1 = nodes(jj(j-1),1);
        lat1 = nodes(jj(j-1),2);
        lon2 = nodes(jj(j),1);
        lat2 = nodes(jj(j),2);

        fprintf(fid,'%8.3f %8.3f \n',[lon1 lat1]);
    end;
    lon1 = nodes(jj(end),1);
    lat1 = nodes(jj(end),2);
    lon2 = nodes(jj(1),1);
    lat2 = nodes(jj(1),2);
    fprintf(fid,'%8.3f %8.3f \n>\n',[lon1 lat1]);

end;

fclose(fid);

namefile = [gmtfile '.names'];
fid=fopen(namefile,'w');


for i=1:M
    eval(['jj=blocks.' char(blocknames(i)) ';']);
    px=mean(nodes(jj,1));
    py=mean(nodes(jj,2));

    fprintf(fid,'%s \n',[num2str(px) ' ' num2str(py) ' 7 0 1 1 ' char(blocknames(i))]);
    
end

fclose(fid);

