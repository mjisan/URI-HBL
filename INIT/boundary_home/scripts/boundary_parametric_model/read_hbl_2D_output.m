function [uv] = read_hbl_2D_output(fi, time, level)
hbl_time = time + 5;
if strcmp(level,'bot')
    u = ncread(fi, 'um_bot', [1,1,hbl_time], [Inf,Inf,1]);
    v = ncread(fi, 'vm_bot', [1,1,hbl_time], [Inf,Inf,1]);
elseif strcmp(level,'top')
    u = ncread(fi, 'um_top', [1,1,hbl_time], [Inf,Inf,1]);
    v = ncread(fi, 'vm_top', [1,1,hbl_time], [Inf,Inf,1]);
end
    
uv = sqrt(u.^2 + v.^2);
end

% [first one] is where to start reading the variable [is how many places to
% read]
