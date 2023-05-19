function [fig_title] = make_plot_rmw(uv, lat, long, plot_title, rmw, mw, xmin, xmax, ymin, ymax, time, x_center, y_center,fig_num, map_on)

mw_txt = sprintf('Max Wind: %0.2f (m/s) Radius of Max Wind: %0.2f (km)', mw, rmw);
%rmw_txt = sprintf('Radius of Mean Wind: %0.2f (km)', rmw);
%info = {mw_txt rmw_txt};
info ={mw_txt};

% Initilize Plot:
fig_title = figure(fig_num);
clf(fig_title)

hold on

%Plot the wind
pcolor(long, lat, uv'); shading flat
colormap jet

%Plot the map
if map_on == 1
       %read coast line
       load('na1.mat');
       load('EastCoast_shoreline.mat')
       coast_lon = EastCoast.lon;
       coast_lat = EastCoast.lat;

    %Plot the coast
    coast_line = line(na1(:,1), na1(:,2), 'color', 'k');
    set([coast_line],'LineWidth',2)

    coast_line = line(coast_lon, coast_lat, 'color', 'k');
    set([coast_line],'LineWidth',2)

    %Plot the track
    line(x_center(1:time), y_center(1:time), 'color', 'k')
    plot(x_center(time),y_center(time),'ro', 'markers', 12)
    line(coast_lon, coast_lat, 'color', 'k')
end

% Add radius of max wind markers
%plot(long(171),lat(162),'ko', 'markers', 12)
%plot(long(181),lat(162),'ko', 'markers', 12)
%plot(long(191),lat(162),'ko', 'markers', 12)

cb = colorbar;
%caxis([0 30]) %incase you want to set it to something specific

%Format the plot titles and stuff
set(gca,'fontsize',18)
title(plot_title)
ylabel(cb, 'Wind Speed (m/s) at 10m')
xlabel('Longitude')
ylabel('Latitude')

% xlim([long(2), long(321)])
% ylim([lat(2), lat(321)])

xlim([long(xmin),long(xmax)])
ylim([lat(ymin),lat(ymax)])

%xlim([-73.54, -70.2218])
%ylim([28.5611, 31.4389])

%xlim([-75,-69])
%ylim([39,43])

%To look at naragansett bay specificaly 
%xlim([-71.8,-71])
%ylim([41,41.9])

%To look at mystic specificaly 
%xlim([-72.1,-71.8])
%ylim([41,41.45])


%.92 for normal
txt = annotation('textbox', [.15, .92, 1, 0],'string', info);
box = txt.LineStyle;
txt.LineStyle = 'none';
txt.FontSize = 12;
txt.Color = 'w';

colormap(jet) 
hold off
end

