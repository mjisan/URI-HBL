function [u,v] = read_hbl_1d_profile(fi, time, i, j)
u = ncread(fi, 'um', [i,j,1,time], [1,1,Inf,1]);
v = ncread(fi, 'vm', [i,j,1,time], [1,1,Inf,1]);
end

% [first one] is where to start reading the variable [is how many places to
% read]
