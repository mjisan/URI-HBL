function [uz,vz,wz,turz] = read_hbl_2d_profile(fi, time,grid_half)

uz = ncread(fi, 'um', [grid_half,grid_half,1,time], [Inf,1,Inf,1]);
vz = ncread(fi, 'vm', [grid_half,grid_half,1,time], [Inf,1,Inf,1]);
wz = ncread(fi, 'wm', [grid_half,grid_half,1,time], [Inf,1,Inf,1]);
turz = ncread(fi, 'tur', [grid_half,grid_half,1,time], [Inf,1,Inf,1]);
end

