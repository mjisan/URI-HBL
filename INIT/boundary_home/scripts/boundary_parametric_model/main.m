%%% Main Plot Storm Script
% Catrina Nowakowski
% 08-18-2018
clear;clc;

%%% Requires these functions to run:
% read_hbl_2D_output.m,read_par.m, read_details.m read_hbl_2d_profile.m
% calc_mw.m, lat_long.m make_plot_rmw.m read_hbl_1d_profile.m
% EastCoast_shoreline.mat na1.mat boundary_parametric.nc boundary_model.nc

%%% Requires diag_table to have:
%um, vm, wm, tur, bot,top...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logic to turn on and off plots and calculations

%%%%%%%%%%%%%%%%%%%%%%%
% Kind of plots to make
parametric_bot = 0; 
parametric_top = 0;
model_bot = 0;
model_top = 1;
vertical_profiles_2d = 1;
vertical_profiles_1d = 1;

% Model resolution (1 or 4 or 10 or 500)
res = 10;
n_rmw = 20; % the domaine you will see the plot on 
nz = 101;
storm_name = 'Carol';
time = 50; %Bob makes land fall at 167 !Rhody land fall 43 !hour after land fall Rhody 47, Irma 80
% For time rember that read_hbl_2D_output adds 5 to it because of the
% spinup steps in the begining of the file, MAKE SURE TO CHANGE it if you
% change the spinup time...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of file locations 
% fi_par = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/CAROL/4KM/OUTPUT/boundary_parametric.nc';
% fi_hbl = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/CAROL/4KM/OUTPUT/boundary_model.nc';

fi_par = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/IDEAL/10KM/OUTPUT/boundary_parametric.nc';
fi_hbl = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/IDEAL/10KM/OUTPUT/boundary_model.nc';

% fi_par = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/CAROL/10KM/OUTPUT_08-30/boundary_parametric.nc';
% fi_hbl = 'http://tds.renci.org:8080/thredds/dodsC/dhs-crc-unc/catrinan/CAROL/10KM/OUTPUT_08-30/boundary_model.nc';

%%%%%%%%%%%%%%%%%%%%%%% Finished setting up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different Cases:
if strcmp(storm_name,'Rhody')
st_time = datetime([2020,08,18,00,00,00]);
map_on = 1;
elseif strcmp(storm_name,'Ideal')
st_time = datetime([1995,01,05,00,00,00]); 
map_on = 0;
elseif strcmp(storm_name,'Storm')
st_time = datetime([1995,01,05,00,00,00]); 
map_on = 0;
elseif strcmp(storm_name,'Irma')
st_time = datetime([2017,09,07,00,00,00]); 
map_on = 1;
elseif strcmp(storm_name,'Bob')
st_time = datetime([1991,09,18,00,00,00]);
map_on = 1;
elseif strcmp(storm_name,'Carol')
st_time = datetime([1954,08,31,00,00,00]);
map_on = 1;
else
    sprintf('Please pick a different storm_name!')
    return
end

if parametric_bot == 1
    model_type = '10m Wind Parametric';
    level = 'bot';
elseif parametric_top == 1
    model_type = '3km Wind Parametric';
    level = 'top';
elseif model_bot == 1
    model_type = '10m Wind Model';
    level = 'bot'; 
elseif model_top == 1
    model_type = '3km Wind Model';
    level = 'top';
end

%model resolution
stepz = 30e-3;
if res == 4
stepx = 4e3;
stepy = 4e3;
res = 1/24;
grid_half = 161;

x_1RMW = 171;
y_1RMW = 161;

x_2RMW = 181;
y_2RMW = 161;

x_3RMW = 191;
y_3RMW = 161;

elseif res == 1
stepx = 1e3;
stepy = 1e3;
elseif res == 10
stepx = 10e3;
stepy = 10e3;

grid_half = 41;

x_1RMW = 45;
y_1RMW = 41;
     
x_2RMW = 49;
y_2RMW = 41;
    
x_3RMW = 53;
y_3RMW = 41;
elseif res == 500
stepx = 500;
stepy = 500;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format Time
time_step = 15; %minutes
time_fmt_single = dateshift(st_time, 'start', 'minute', time*time_step);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autogenerate file name:
file_name = sprintf('%s %s.png',storm_name, model_type );
plot_title = sprintf('%s %s %s',storm_name, model_type, datestr(time_fmt_single) );


fig_num = 1; %Dont change this, it's a counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the data for what ever is turned on

%%%%%%%%%%%%%%%%%%%%%%%%%%
% read parameters from par
fprintf('reading parameters...')
fi = fi_par;
[x,y,x_center, y_center] = read_details(fi);
%format lat and lon
[lat, long ] = lat_long( y_center(time),x_center(time), ...
                        length(y), length(x), stepy, stepx);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% read for specific plots
if parametric_bot == 1 
    fi = fi_par;
    [uv,mask] = read_par(fi, time, level);
    uv = uv;
end
if parametric_top == 1
    fi = fi_par;
    [uv,mask] = read_par(fi, time, level);
    uv= uv;
end
if model_bot == 1
    fi = fi_hbl;
    [uv] = read_hbl_2D_output(fi, time, level);
    uv = uv/1.1324; %as of 4-15-2018
end  
if model_top == 1
    fi = fi_hbl;
    [uv] = read_hbl_2D_output(fi, time, level);
    uv = uv/1.1324; %as of 4-15-2018
end
if vertical_profiles_2d == 1
    fi = fi_hbl;
    [uz,vz,wz,turz] = read_hbl_2d_profile(fi, time,grid_half);
end
if vertical_profiles_1d == 1
    fi = fi_hbl;

    [u,v] = read_hbl_1d_profile(fi, time, x_1RMW, y_1RMW);
    um_1RMW(1,:) = u(1,1,:);
    vm_1RMW(1,:) = v(1,1,:);

    [u,v] = read_hbl_1d_profile(fi, time, x_2RMW, y_2RMW);
    um_2RMW(1,:) = u(1,1,:);
    vm_2RMW(1,:) = v(1,1,:);

    [u,v] = read_hbl_1d_profile(fi, time, x_3RMW, y_3RMW);
    um_3RMW(1,:) = u(1,1,:);
    vm_3RMW(1,:) = v(1,1,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Horizontal Plot

if parametric_bot == 1 || parametric_top == 1 || model_bot == 1 || model_top == 1
    %calc the max wind and RMW
    [mw, rmw, xmin, xmax, ymin, ymax] = calc_mw(x,y,stepx,stepy, uv, n_rmw);    
    % Make plot at that step
    fig_title = make_plot_rmw(uv, lat, long, plot_title, rmw, mw, xmin, xmax, ymin, ymax, time, x_center, y_center,fig_num, map_on);
    fig_num = fig_num + 1;
    % Print the plot
    print(fig_title, file_name, '-dpng')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if vertical_profiles_2d == 1
    grid_indx = 42 ; %162; %42
    
    z_axis = (0:1:(nz-1))*stepz;
    x_axis = (0:1:(grid_indx-1))*stepx;
    y_axis = (0:1:(grid_indx-1))*stepy;
    
    
    %Define colormap
    c1=[0 0 1]; %G
    c2=[0 1 1]; %c
    c3=[0 1 0]; %b
    n1=20;
    n2=20;
    cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];

    umz(:,:) = uz(:,1,:);
    uz_levels = -2:-2:-10;
    
    figure(fig_num)
    fig_num = fig_num + 1;
    subplot(4,1,1)
    hold on 
    [x,z] = meshgrid((x_axis(1,1:grid_indx)/40e3),z_axis);
    pcolor( x,z,umz'); shading flat
    contour(x,z,umz',uz_levels,'k')
    xlim([0 5])
    cb = colorbar;
    colormap(cmap')
    caxis([-10 2])
    xlabel('Radial Wind (m/s)')
    ylabel('Z (km)') 
    hold off

    %Define colormap
    c1=[0 1 0]; %G
    c2=[1 1 0]; %Y
    c3=[1 0 0]; %R
    n1=20;
    n2=20;
    cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];

    vmz(:,:) = vz(:,1,:);
    vz_levels = 0:10:50;
    
    subplot(4,1,2)
    hold on 
    [y,z] = meshgrid(y_axis(1,1:grid_indx)/40e3,z_axis);
    pcolor( y,z,vmz'); shading flat
    contour(y,z,vmz',vz_levels,'k')
    xlim([0 5])
    cb = colorbar;
    colormap(cmap')
    caxis([0 50])
    xlabel('Azimuthal Wind (m/s)')
    ylabel('Z (km)') 
    hold off

    wmz(:,:) = wz(:,1,:);
    wz_levels = 0:0.05:0.2;
    
    subplot(4,1,3)
    hold on 
    [y,z] = meshgrid(y_axis(1,1:grid_indx)/40e3,z_axis);
    pcolor( y,z,wmz'); shading flat
    contour(y,z,wmz',wz_levels,'k')
    xlim([0 5])
    cb = colorbar;
    colormap(cmap')
    caxis([0 0.2])
    xlabel('Vertical Wind (m/s)')
    ylabel('Z (km)') 
    hold off

    turmz(:,:) = turz(:,1,:);
    turz_levels = 0:20:60;
    
    subplot(4,1,4)
    hold on 
    [y,z] = meshgrid(y_axis(1,1:grid_indx)/40e3,z_axis);
    pcolor( y,z,turmz'); shading interp
    contour(y,z,turmz',turz_levels,'k')
    xlim([0 5])
    cb = colorbar;
    colormap(cmap')
    caxis([0 60])
    xlabel('Turbulent Viscosity(m2/s)')
    ylabel('Z (km)') 
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if vertical_profiles_1d == 1
    z_axis = (0:1:(nz-1))*stepz;
    x_axis = (0:1:(length(x)))*stepx;
    y_axis = (0:1:(length(y)))*stepy;


    figure(fig_num)
    subplot(1,2,1)
    fig_num = fig_num + 1;
    hold on
    plot(um_1RMW,z_axis, '-')
    plot(um_2RMW,z_axis, '--')
    plot(um_3RMW,z_axis, '-.')
    xlabel('Radial Wind (m/s)')
    ylabel('z (km)')
    xlim([-10 5])
    hold off

    subplot(1,2,2)
    hold on
    plot(vm_1RMW,z_axis, '-')
    plot(vm_2RMW,z_axis, '--')
    plot(vm_3RMW,z_axis, '-.')
    xlabel('Azimuthal Wind (m/s)')
    ylabel('z (km)')
    xlim([10 50])
    hold off
end








%%% IF you want to take out low values from the plot

    %    for i = 1:length(uv(:,1))
    %    for j = 1:length(uv(1,:))
    %        
    %        if uv(i,j) < 10
    %            uv(i,j) = nan;
    %        end
    %        
    %    end
    %    end