function fnums=GetFaultsConnected(faults,ifnum)
% fnums=GetFaultsConnected(faults,ifnum)

nset = unique(faults(ifnum,1:2));
growing=1;
fnums=ifnum;
% it=0;
while growing
%     it=it+1;
%     disp(num2str(it));
    
    fcon1=[];
    fcon2=[];
    for k=1:length(nset)
        fn1=find(nset(k)==faults(:,1));
        fn2=find(nset(k)==faults(:,2));
        fcon1=unique([fcon1;fn1]);
        fcon2=unique([fcon2;fn2]);
    end
    
    fnew = setdiff(unique([fcon1;fcon2]),fnums);
    if ~isempty(fnew)
       fnums=[fnums;fnew];
       nset = unique(faults(fnums,1:2));
    else
        growing=0;
    end
    
%     disp(fnums)
%     pause
end





