function [mw, rmw, xmin, xmax, ymin, ymax] = calc_mw(x,y,stepx,stepy, uv, n_rmw)

%Set the constants
nx=length(x);ny=length(y);

%This makes the grid for when I want the center == 0
%the +1 in the first parts are to adjust for the fact that the grid center
%is at (151,151)
xgrid_center = (-floor(nx/2)):1:(floor(nx/2)-1);
xgrid_center = xgrid_center *  stepx * 1.e-3;

ygrid_center = (-floor(ny/2)):1:(floor(ny/2)-1);
ygrid_center = ygrid_center *  stepy * 1.e-3;

mw = max(max(uv(:,:)));
[mw_x,mw_y] = find(uv(:,:) == mw);
rmw = sqrt((xgrid_center(mw_x(1))^2) + (ygrid_center(mw_y(1))^2));
% rmw = 40;
n_rmw = n_rmw*rmw;

[gar, xmin] = min(abs(xgrid_center + n_rmw));
[gar, xmax] = min(abs(xgrid_center - n_rmw));

[gar, ymin] = min(abs(ygrid_center + n_rmw));
[gar, ymax] = min(abs(ygrid_center - n_rmw));
end

