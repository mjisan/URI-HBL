function [x,y,x_center, y_center] = read_details(fi)
%Read parameters from par
% x,y,x_center, y_center
x = ncread(fi, 'x_axis');
y = ncread(fi, 'y_axis');
x_center = ncread(fi, 'x_center');
y_center = ncread(fi, 'y_center');
end

