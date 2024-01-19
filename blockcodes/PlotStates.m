function PlotStates(opt,opt2)
%
%  plots outlies of western states CA/NV/OR/ID
%
% if opt==1 then use 0<lon<360
% else use lon -180<lon<180
%
% if opt2=1 plots western states
% if opt2= other plots lower 48 states

dlon=0;
% if nargin==1
    if opt==1
        dlon = 360;
    end
% end;

states = shaperead('usastatehi');
if (opt2==1)
   iuse = [3 5 6 12 28 31 37 26 47 44 50];
else
   iuse = 1:51;
end
for i=iuse
    x = getfield(states(i),'X');
    y = getfield(states(i),'Y');
    plot(x+dlon,y,'k-');
    hold on;
end


    
    
