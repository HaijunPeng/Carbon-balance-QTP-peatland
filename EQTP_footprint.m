%Hongyuan Peatland EC footprint analysis
%Jan 21, 2021
%J. Chi

clear all; clc; close all;
warning off;
addpath 'C:\Users\jihi0001\Documents\MATLAB';
mydat = readtable('C:/Users/jihi0001/Dropbox/Publications/Peng_Chi_Manuscript/Data/GHG_EQTP_lv1.csv');
z = 2.5; %EC height above the ground, m
hc = 0.1; %canopy height, 10 cm, constant throughout the year
d = 0.67*hc;
zm = z-d;
z0 = 0.1*hc; % roughness length, 10% of the canopy height
r = 10:10:90;
umean=NaN;

ol_day = mydat.L_m(mydat.Rg_Wm2 > 10); %daytime
ol_night = mydat.L_m(mydat.Rg_Wm2 < 10); %nighttime
mydat.sigmav = sqrt(mydat.v_var_m2s2);
sigmav_day = mydat.sigmav(mydat.Rg_Wm2 > 10);
sigmav_night = mydat.sigmav(mydat.Rg_Wm2 < 10);
ustar_day = mydat.ustar_ms(mydat.Rg_Wm2 > 10);
ustar_night = mydat.ustar_ms(mydat.Rg_Wm2 < 10);
wind_dir_day = mydat.wind_dir_deg(mydat.Rg_Wm2 > 10);
wind_dir_night = mydat.wind_dir_deg(mydat.Rg_Wm2 < 10);
h_day = repmat(800, size(ol_day)); %daytime PBL 800m averaged for southeast QTP
h_night = repmat(200, size(ol_night)); %nighttime PBL 200m averaged for southeast QTP

[FFP1, flag_error1] = calc_footprint_FFP_climatology(zm, z0, umean, h_day, ol_day ,sigmav_day, ustar_day, wind_dir_day,...
    'domain', [-300 300 -300 300], 'dx', 5, 'dy', 5, 'r', r);

[FFP2, flag_error2] = calc_footprint_FFP_climatology(zm, z0, umean, h_night, ol_night ,sigmav_night, ustar_night, wind_dir_night,...
    'domain', [-300 300 -300 300], 'dx', 5, 'dy', 5, 'r', r);

figure(1)
surf(FFP1(1).x_2d, FFP1(1).y_2d, FFP1(1).fclim_2d); shading flat; view(2);
hold all;
for i=1:length(r)
z = FFP1(i).fr.*10.*ones(size(FFP1(i).yr));
plot3(FFP1(i).xr, FFP1(i).yr, z, 'r')
end

figure(2)
surf(FFP2(1).x_2d, FFP2(1).y_2d, FFP2(1).fclim_2d); shading flat; view(2);
hold all;
for i=1:length(r)
z = FFP2(i).fr.*10.*ones(size(FFP2(i).yr));
plot3(FFP2(i).xr, FFP2(i).yr, z, 'r')
end


R = 6371007.181;  % the radius of the idealized sphere representing the Earth, m
lat = 32.761187;
lon = 102.543988; 

dn1 = FFP1(3).yr';
de1 = FFP1(3).xr';
% Coordinate offsets in radians
dLat1 = dn1./R;
dLon1 = de1./(R*cos(pi*lat/180));

% % OffsetPosition, decimal degrees
latO1 = lat + dLat1*180./pi;
lonO1 = lon + dLon1*180./pi;

coord1 = [FFP1(3).xr', FFP1(3).yr', lonO1, latO1];
fid1 = fopen('EQTP_DaytimeFootprint_30p.csv', 'w'); %update
fprintf(fid1,'%s\n', 'x, y, long, lat');
fclose(fid1);
% %write data to end of file
dlmwrite('EQTP_DaytimeFootprint_30p.csv', coord1, '-append', 'precision', '%7.5f'); %update

dn2 = FFP2(3).yr';
de2 = FFP2(3).xr';
% Coordinate offsets in radians
dLat2 = dn2./R;
dLon2 = de2./(R*cos(pi*lat/180));

% % OffsetPosition, decimal degrees
latO2 = lat + dLat2*180./pi;
lonO2 = lon + dLon2*180./pi;

coord2 = [FFP2(3).xr', FFP2(3).yr', lonO2, latO2];
fid2 = fopen('EQTP_NighttimeFootprint_30p.csv', 'w'); %update
fprintf(fid2, '%s\n', 'x, y, long, lat');
fclose(fid2);
% %write data to end of file
dlmwrite('EQTP_NighttimeFootprint_30p.csv', coord2, '-append', 'precision', '%7.5f'); %update