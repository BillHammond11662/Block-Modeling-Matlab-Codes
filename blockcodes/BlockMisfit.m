function [chi2dof,dof,n,rmsn,rmse,rmstot]=BlockMisfit(blocks,nodes,lon,lat,ve,vn,se,sn,vepred,vnpred,flag)
% function [chi2dof,dof,n,rmsn,rmse,rmstot]=BlockMisfit(blocks,nodes,lon,lat,ve,vn,se,sn,vepred,vnpred,flag)
%
% if flag==1 then send text to screen

blocknames=fieldnames(blocks);
M=length(blocknames);

if flag==1
    disp('          Block    RMSE    RMSN  Chi^2/dof  NumSites');
end

modelnodes=struct2cell(blocks);
chi2dof=nan(M,1);

dof = nan(M,1);
rmse=nan(M,1);
rmsn=nan(M,1);
rmstot=nan(M,1);
n=nan(M,1);

% keyboard;

for i=1:M
    
    j=modelnodes{i};
    k = find(inpolygon(lon,lat,nodes(j,1),nodes(j,2)));
    n(i) = length(k);
    
    dof(i) = 2*length(k) - 3;
    
    rmse(i) = rms(ve(k) - vepred(k));
    rmsn(i) = rms(vn(k) - vnpred(k));
    rmstot(i) = rms([vn(k) - vnpred(k);vn(k) - vnpred(k)]);

    chi2 = sum(((vepred(k)-ve(k))./se(k)).^2,'omitnan') + sum(((vnpred(k)-vn(k))./sn(k)).^2,'omitnan');

    if dof(i)>0    
        chi2dof(i,1)=chi2/dof(i);
        str = [sprintf('%15s',char(blocknames(i))) '  ' sprintf('%6.2f',rmse(i)) '  ' sprintf('%6.2f',rmsn(i)) '  ' sprintf('%8.1f',chi2dof(i))  ' ' sprintf('%6.0f',n(i)) ];
    else
        dof(i)=NaN;
        str = [sprintf('%15s',char(blocknames(i))) '   NaN    NaN   NaN     ' sprintf('%3.0f',n(i)) ];
    end

    if flag==1
        disp(str);
    end

%     if isnan(chi2dof(i))
%         keyboard;
%     end
end

% keyboard;
