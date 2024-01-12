
figure(8);
clf

rese=ve(ilist)-vepred(ilist)*1000;  % in mm/yr
resn=vn(ilist)-vnpred(ilist)*1000;

nrese=rese./se(ilist);
nresn=resn./sn(ilist);

vemad = median(abs(rese-median(rese,'omitnan')),'omitnan');
vnmad = median(abs(resn-median(resn,'omitnan')),'omitnan');

vemadnorm = median(abs(nrese-median(nrese,'omitnan')),'omitnan');
vnmadnorm = median(abs(nresn-median(nresn,'omitnan')),'omitnan');


subplot(223);
histogram(nrese)
title({'East Residuals Normalized by Uncertainty',['NRMSE = ' num2str(rms(nrese))],['MADE  = ' num2str(vemadnorm)]});

subplot(224);
histogram(nresn);
title({'North Residuals Normalized by Uncertainty',['NRMSN = ' num2str(rms(nresn))],['MADN = ' num2str(vnmadnorm)]});

subplot(221);
histogram(rese)
xlabel('mm/yr');
title({'East Residuals',['RMSE = ' num2str(rms(rese))],['MADE  = ' num2str(vemad)]});

subplot(222);
histogram(resn);
xlabel('mm/yr');
title({'North Residuals',['RMSN = ' num2str(rms(resn))],['MADN = ' num2str(vnmad)]});

orient portrait

rmsfile = 'nrms.txt';

fid=fopen(rmsfile,'w+');
fprintf(fid,'rese=%6.4f\n',rms(rese));
fprintf(fid,'resn=%6.4f\n',rms(resn));
fprintf(fid,'nrese=%6.4f\n',rms(nrese));
fprintf(fid,'nresn=%6.4f\n',rms(nresn));
fprintf(fid,'nrms=%6.4f\n',rms([nrese;nresn]));
fprintf(fid,'rms=%6.4f\n',rms([rese;resn]));

fprintf(fid,'vemad=%6.4f\n',vemad);
fprintf(fid,'vnmad=%6.4f\n',vnmad);
fprintf(fid,'vemadnorm=%6.4f\n',vemadnorm);
fprintf(fid,'vnmadnorm=%6.4f\n',vnmadnorm);

fprintf(fid,'chi2=%6.4f\n',chi2);
fprintf(fid,'dof=%6.4f\n',dof);
fprintf(fid,'chi2dof=%6.4f\n',chi2/dof);
fclose(fid);