function c = convexx(X)
% Array c is 1 where polygon is convex, -1 where curve is concave
%  Input: X = n by 2 array of x and y values
%
% Modified from:
% Are Mjaavatten (2024). Identify convex and concave parts of polygon
% (https://www.mathworks.com/matlabcentral/fileexchange/102169-identify-convex-and-concave-parts-of-polygon), 
% MATLAB Central File Exchange. Retrieved January 16, 2024. 
% Are Mjaavatten, November 2021
% 

flag=0;
if (X(end,1)==X(1,1) && X(end,2)==X(1,2))
    X(end,:)=[];
    flag=1;
end

n = size(X,1);
Z = [X,zeros(n,1)];
c = NaN(n,1);
for i = 1:n
    if i==1
        cc = cross(Z(i,:)-Z(end,:),Z(i+1,:)-Z(end,:));
        c(i) = -sign(cc(3));
    elseif i==n
        cc = cross(Z(i,:)-Z(i-1,:),Z(1,:)-Z(i-1,:));
        c(i) = -sign(cc(3));
    else
        cc = cross(Z(i,:)-Z(i-1,:),Z(i+1,:)-Z(i-1,:));
        c(i) = -sign(cc(3));
    end
end

% Change sign if the points are not in clockwise order:
clockwise = sign(sum(circshift(X(:,1),-1).*X(:,2)-X(:,1).*circshift(X(:,2),-1)));
c = c*clockwise;

if flag==1
    c=[c;c(1)];
end