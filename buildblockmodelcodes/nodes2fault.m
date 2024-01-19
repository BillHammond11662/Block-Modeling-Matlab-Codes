function fnum = nodes2fault(faults,n1,n2)
% fnum = nodes2fault(faults,n1,n2)

fnum= find((faults(:,1)==n1 & faults(:,2)==n2) | ...
           (faults(:,1)==n2 & faults(:,2)==n1));

    


