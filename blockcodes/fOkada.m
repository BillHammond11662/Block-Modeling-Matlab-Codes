function [f]=fOkada(xi,eta,z,q,dip,alpha,typ,term,i)

R = sqrt(xi.^2 + eta.^2 + q.^2);
ytild = eta*cos(dip) + q*sin(dip);
dtild = eta*sin(dip) - q*cos(dip);
ctild = dtild + z;
if abs(q)>eps
   %theta = atan2(xi*eta,q*R);
   theta = atan(xi*eta/(q*R));
else
   theta =0;
end

X=sqrt(xi.^2 + q.^2);
 
if abs(R+xi)>eps
    X11 = 1/(R.*(R+xi));
    X32 = (2*R + xi)./(R.^3*((R+xi).^2));
    X53 = (8*R.^2 + 9*R.*xi + 3*xi.^2)./((R.^5).*((R+xi).^2));
else
    X11=0;
    X32=0;
    X53=0;
end

% figure(2);
% plot(xi,X11,'go');
% hold on;

if abs(R+eta)>eps
    Y11 = 1./(R.*(R+eta));
    Y32 = (2*R + eta)./(R.^3*((R+eta).^2));
    Y53 = (8*R.^2 + 9*R.*eta + 3*eta.^2)./(R.^5 .* ((R+eta).^2) ) ;
else
    Y11=0;
    Y32=0;
    Y53=0;
end

h = q*cos(dip) - z;
Z32 = sin(dip)./R.^3 - h.*Y32;
Z53 = 3*sin(dip)./R.^5 - h.*Y53;

Y0 = Y11 - xi.^2.*Y32;
Z0 = Z32 - xi.^2.*Z53;

if abs(cos(dip))>eps
   if abs(R+eta)>eps
      I3 = (1/cos(dip))*(ytild./(R + dtild)) - (1/cos(dip).^2)*(log(R+eta) - sin(dip)*log(R+dtild));
   else
      I3 = (1/cos(dip))*(ytild./(R + dtild)) + (1/cos(dip).^2)*(log(R-eta) - sin(dip)*log(R+dtild));
   end
   I4 = (sin(dip)/cos(dip))*(xi./(R + dtild)) + ...
       (2/(cos(dip).^2))*atan((eta.*(X+q*cos(dip)) + X.*(R+X)*sin(dip))./(xi.*(R+X)*cos(dip)));
else
   if abs(R+eta)>eps
      I3 = .5*(eta./(R + dtild) + ytild.*q./((R + dtild).^2) - log(R+eta));
   else
      I3 = .5*(eta./(R + dtild) + ytild.*q./((R + dtild).^2) + log(R-eta));
   end
   I4 = .5*(xi*ytild./(R+dtild).^2);
end
if abs(xi)<eps
    I4=0;
end

I1 = (-xi./(R + dtild))*cos(dip) - I4*sin(dip);
I2 = log(R + dtild) + I3*sin(dip);


if strcmp(typ,'strike')

    if term=='A'
        if i==1
            f = theta/2 + (alpha/2).*xi.*q.*Y11;
        elseif i==2
            f = (alpha/2).*(q./R);
        elseif i==3
            if abs(R+eta)>eps
               f = ((1-alpha)/2)*log(R+eta) - (alpha/2)*(q.^2)*Y11;
            else
               f =-((1-alpha)/2)*log(R-eta) - (alpha/2)*(q.^2)*Y11;
            end
        end

    elseif term=='B'
        if i==1
            f = -xi*q*Y11 - theta - ((1-alpha)/alpha)*I1*sin(dip);
        elseif i==2
            f = -q/R               + ((1-alpha)/alpha)*(ytild/(R+dtild))*sin(dip);
        elseif i==3
            f = (q.^2)*Y11         - ((1-alpha)/alpha)*I2*sin(dip);
        end

    elseif term=='C'
        if i==1
            f = (1-alpha)*xi*Y11*cos(dip)                 - alpha*xi*q*Z32;
        elseif i==2
            f = (1-alpha)*(cos(dip)/R + 2*q*Y11*sin(dip)) - alpha*ctild*q/(R.^3);
        elseif i==3
            f = (1-alpha)*q*Y11*cos(dip)                  - alpha*(ctild*eta/(R.^3) - z*Y11 + xi.^2*Z32);
        end
    end

end

if strcmp(typ,'dip')

    if term=='A'
        if i==1
            f = (alpha/2)*(q/R);
        elseif i==2
            f = theta/2 + (alpha/2)*eta*q*X11;
        elseif i==3
            if abs(R + xi)>eps
               f = ((1-alpha)/2)*log(R+xi) - (alpha/2)*(q.^2)*X11;
            else
               f =-((1-alpha)/2)*log(R-xi) - (alpha/2)*(q.^2)*X11;
            end
        end

    elseif term=='B'
        if i==1
            f = -q/R               + ((1-alpha)/alpha)*I3*sin(dip)*cos(dip);
        elseif i==2
            f = -eta*q*X11 - theta - ((1-alpha)/alpha)*(xi/(R+dtild))*sin(dip)*cos(dip);
%             figure(2);
%             subplot(131);plot3(xi,eta,-eta*q*X11,'ro');
%             hold on;
%             subplot(132);plot3(xi,eta,-theta,'go'); hold on;
%             subplot(133);plot3(xi,eta,-((1-alpha)/alpha)*(xi/(R+dtild))*sin(dip)*cos(dip),'bo'); hold on;
        elseif i==3
            f = (q.^2)*X11        + ((1-alpha)/alpha)*I4*sin(dip)*cos(dip);
        end

    elseif term=='C'
        if i==1
            f = (1-alpha)*(cos(dip)/R) - q*Y11*sin(dip)  - alpha*ctild*q/(R.^3);
        elseif i==2
            f = (1-alpha)*ytild*X11                      - alpha*ctild*eta*q*X32;
        elseif i==3
            f = -dtild*X11 - xi*Y11*sin(dip)             - alpha*ctild*(X11 - q.^2*X32);
        end
    end

end


if strcmp(typ,'tensile')

    if term=='A'
        if i==1
            if abs(R+eta)>eps
               f = -((1-alpha)/2)*log(R+eta) - (alpha/2)*(q.^2)*Y11;
            else
               f =  ((1-alpha)/2)*log(R-eta) - (alpha/2)*(q.^2)*Y11;
            end
        elseif i==2
            if abs(R + xi)>eps
               f = -((1-alpha)/2)*log(R+xi)  - (alpha/2)*(q.^2)*X11;
            else
               f =  ((1-alpha)/2)*log(R-xi)  - (alpha/2)*(q.^2)*X11;
            end
        elseif i==3
            f = theta/2 - (alpha/2)*q*(eta*X11 + xi*Y11);
        end

    elseif term=='B'
        if i==1
            f = Y11*q.^2             - ((1-alpha)/alpha)*I3*(sin(dip).^2);
        elseif i==2
            f = X11*q.^2             + ((1-alpha)/alpha)*(xi/(R+dtild))*(sin(dip).^2);
        elseif i==3
            f = q*(eta*X11 + xi*Y11) - theta - ((1-alpha)/alpha)*I4*(sin(dip).^2);
        end

    elseif term=='C'
        if i==1
           f = -(1-alpha)*(sin(dip)/R + q*Y11*cos(dip)) - alpha*(z*Y11 - Z32*q.^2);
        elseif i==2
           f =  (1-alpha)*2*xi*Y11*sin(dip) + dtild*X11 - alpha*ctild*(X11 - X32*q.^2);
        elseif i==3
           f =  (1-alpha)*(ytild*X11 + xi*Y11*cos(dip)) + alpha*q*(ctild*eta*X32 - Z32*xi);
        end
    end

end