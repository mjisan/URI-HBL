function [uv,mask] = read_par(fi, time, level)
if strcmp(level,'bot')
    u = ncread(fi, 'u10', [1,1,time], [Inf,Inf,1]);
    v = ncread(fi, 'v10', [1,1,time], [Inf,Inf,1]);
elseif strcmp(level,'top')
    u = ncread(fi, 'utop', [1,1,time], [Inf,Inf,1]);
    v = ncread(fi, 'vtop', [1,1,time], [Inf,Inf,1]);
end

mask = ncread(fi, 'land_mask', [1,1,time], [Inf,Inf,1]);

uv = sqrt(u.^2 + v.^2);
end

