
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