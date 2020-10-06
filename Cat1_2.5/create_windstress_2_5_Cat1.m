% This script caculate surface wind speed and air pressure of idealized
% hurricane.The formula can be found in Holland(1980). Wind speed is ramped
% up over one day.
%
% Translation speed = 2.5 m/s, or 9 km/h
% Maximum Wind Speed = 37.54 m/s, or 135.1 km/h (Cat 1)
clc
clear
% set output file names
uwndfn = 'ocean_frc_uwstress_2_5_Cat1.nc';
vwndfn = 'ocean_frc_vwstress_2_5_Cat1.nc';
% load grid informations
load grid_info13.mat %8 km grid
% set wind forcing time (Unit:day). time_start could be the start time of
% your simulation and time_end could be the end time.
time_start = 0;
time_end = 10;
time_interval = 1/24;
wind_time = time_start:time_interval:time_end;
% set the initial locations of hurricane center
x_center0 = 3*1e6;
y_center0 = 1.5*1e6;
% set moving speed of hurricane center (Unit: km hour-1)
center_s_x = -9.0; %2.5 m/s
center_s_y = 0;
%
% P_C is the pressure in hurricane center and R_Max is the radius of
% maximum wind. Default values are from Fig.2 of Holland(1980)
P_N = 1005; 
P_C = 975; % Cat 1 is > 979
R_Max = 50*1e3; %unit:m
% cff_B is the parameter B in ther formula
cff_B = 1.5 + (980 - P_C) / 120;
% Air density (Unit: kg m-3)
rho_A = 1.15;

% caculate wind speed and air pressure
for tindx = 1:length(wind_time)
   x_center = x_center0 + center_s_x * (wind_time(tindx) - wind_time(1)) * 24 * 1000;
   y_center = y_center0 + center_s_y * (wind_time(tindx) - wind_time(1)) * 24 * 1000;
   center(tindx,1) = x_center;
   center(tindx,2) = y_center;
   if tindx <= 24
       a(tindx) = (tindx-1)/24;
   else
       a(tindx) = 1;
   end
   dist = sqrt((x_rho - x_center) .^ 2 + (y_rho - y_center) .^ 2);
   dist(dist < 1) = 1;
   %
   wind_sp = sqrt(cff_B * (R_Max ./ dist) .^ cff_B .* (P_N - P_C) .* 100 .*  ...
       exp(-(R_Max ./ dist) .^ cff_B) ./ rho_A + dist .^2 .* f .^ 2 / 4) - ...
       dist .* f / 2;
   % 
   %disp(max(max(wind_sp)))
   %
   Cd = zeros(size(wind_sp)) + 1.2/1000;
   Cd(wind_sp > 11 & wind_sp <= 19) = (0.49 + 0.065 * wind_sp(wind_sp > 11 & wind_sp <= 19)) / 1000;
   Cd(wind_sp > 19) = (1.364 + 0.0234 * wind_sp(wind_sp > 19) - 0.0002 * wind_sp(wind_sp > 19) .^ 2) / 1000;
   %
   wind_st = rho_A .* Cd .* wind_sp .^ 2;%for test29, multiply by a(tindx)
   %
   x_tmp = x_rho - x_center;
   y_tmp = y_rho - y_center;
   wind_st_xtmp = wind_st .* (-(y_tmp) ./ dist);
   wind_st_ytmp = wind_st .* (x_tmp ./ dist);
   
   wind_st_x(tindx,:,:) = a(tindx).*interp2(x_rho,y_rho,wind_st_xtmp,x_u,y_u);
   wind_st_y(tindx,:,:) = a(tindx).*interp2(x_rho,y_rho,wind_st_ytmp,x_v,y_v);
   %
   %disp(tindex)
end
pcolor(x_rho,y_rho,wind_sp)
hold on
plot(center(:,1),center(:,2))

% 
% %---------------Write netCDF files-------------------------------------
% %
[Mp,Lp] = size(x_rho);
L = Lp - 1;
M = Mp - 1;
%Create wind stress file
nc_create_empty(uwndfn)
%Add dimensions
nc_add_dimension(uwndfn,'xi_u',L);
nc_add_dimension(uwndfn,'eta_u',Mp);
nc_add_dimension(uwndfn,'xi_v',Lp);
nc_add_dimension(uwndfn,'eta_v',M);
nc_add_dimension(uwndfn,'xi_rho',Lp);
nc_add_dimension(uwndfn,'eta_rho',Mp);
nc_add_dimension(uwndfn,'xi_psi',L);
nc_add_dimension(uwndfn,'eta_psi',M);
nc_add_dimension(uwndfn,'sms_time',length(wind_time));
nc_add_dimension(uwndfn,'x',size(x_rho,2));
nc_add_dimension(uwndfn,'y',size(x_rho,1));
nc_add_dimension(uwndfn,'one',1);
%
%Add variables
var_tmp.Name = 'spherical';
var_tmp.Dimension = {'one'};
nc_addvar(uwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'x';
var_tmp.Dimension = {'eta_u','xi_u'};
nc_addvar(uwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'y';
var_tmp.Dimension = {'eta_u','xi_u'};
nc_addvar(uwndfn,var_tmp)
clear var_tmp
%Add variables
var_tmp.Name = 'sms_time';
var_tmp.Datatype ='double';
var_tmp.Dimension = {'sms_time'};
var_tmp.Attribute(1).Name = 'units';
var_tmp.Attribute(1).Value = 'day';
var_tmp.Attribute(2).Name = 'long_name';
var_tmp.Attribute(2).Value = 'surface momentum stress time';
nc_addvar(uwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'sustr';
var_tmp.Datatype ='double';
var_tmp.Dimension = {'sms_time', 'eta_u', 'xi_u'};
var_tmp.Attribute(1).Name = 'units';
var_tmp.Attribute(1).Value = 'N m-2';
var_tmp.Attribute(2).Name = 'long_name';
var_tmp.Attribute(2).Value = 'surface u-momentum stress';
var_tmp.Attribute(3).Name = 'coordinates';
var_tmp.Attribute(3).Value = 'x y sms_time';
nc_addvar(uwndfn,var_tmp)
%
clear var_tmp
% Write variables
nc_varput(uwndfn,'spherical',0)
nc_varput(uwndfn,'x',x_u);
nc_varput(uwndfn,'y',y_u);
nc_varput(uwndfn,'sms_time',wind_time)
nc_varput(uwndfn,'sustr',wind_st_x)
% %-------------------
%Create V-Wind file
nc_create_empty(vwndfn)
%Add dimensions
nc_add_dimension(vwndfn,'xi_u',L);
nc_add_dimension(vwndfn,'eta_u',Mp);
nc_add_dimension(vwndfn,'xi_v',Lp);
nc_add_dimension(vwndfn,'eta_v',M);
nc_add_dimension(vwndfn,'xi_rho',Lp);
nc_add_dimension(vwndfn,'eta_rho',Mp);
nc_add_dimension(vwndfn,'xi_psi',L);
nc_add_dimension(vwndfn,'eta_psi',M);
nc_add_dimension(vwndfn,'sms_time',length(wind_time));
nc_add_dimension(vwndfn,'lon',size(x_rho,2));
nc_add_dimension(vwndfn,'lat',size(x_rho,1));
nc_add_dimension(vwndfn,'one',1);
%
%Add variables
var_tmp.Name = 'spherical';
var_tmp.Dimension = {'one'};
nc_addvar(vwndfn,var_tmp)
clear var_tmp
%
%
var_tmp.Name = 'x';
var_tmp.Dimension = {'eta_v','xi_v'};
nc_addvar(vwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'y';
var_tmp.Dimension = {'eta_v','xi_v'};
nc_addvar(vwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'sms_time';
var_tmp.Datatype ='double';
var_tmp.Dimension = {'sms_time'};
var_tmp.Attribute(1).Name = 'units';
var_tmp.Attribute(1).Value = 'day';
var_tmp.Attribute(2).Name = 'long_name';
var_tmp.Attribute(2).Value = 'surface momentum stress time';
nc_addvar(vwndfn,var_tmp)
clear var_tmp
%
var_tmp.Name = 'svstr';
var_tmp.Datatype ='double';
var_tmp.Dimension = {'sms_time', 'eta_v', 'xi_v'};
var_tmp.Attribute(1).Name = 'units';
var_tmp.Attribute(1).Value = 'N m-2';
var_tmp.Attribute(2).Name = 'long_name';
var_tmp.Attribute(2).Value = 'surface v-momentum stress';
var_tmp.Attribute(3).Name = 'coordinates';
var_tmp.Attribute(3).Value = 'lon lat sms_time';
nc_addvar(vwndfn,var_tmp)
%
clear var_tmp
% Write variables
nc_varput(vwndfn,'spherical',0)
nc_varput(vwndfn,'x',x_v);
nc_varput(vwndfn,'y',y_v);
nc_varput(vwndfn,'sms_time',wind_time)
nc_varput(vwndfn,'svstr',wind_st_y)
%
disp('Done!')
