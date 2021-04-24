% Martian background atmosphere for CO2 clouds validation: MOLA pass # 260
% Parameters from MCD5.3 are saved in a mat file

% only parameters till 95 km are taken

function[T,rho,p,R,gamma,kvisc,thermdiffus,H,C,U] = Mars_MOLApass260(z_array)

% Ls 316.4deg. Latitude 84.0N Longitude 72.5E Local time 0.0h 
load('atmo_pass260.mat') % this mat file must be provided (created from MCD data)

% read MCD data
z_data = table2array(atmo(:,1)); % m
T_data = table2array(atmo(:,2)); % K
rho_data = table2array(atmo(:,3)); % kg/m3
P_data = table2array(atmo(:,4)); % Pa
U_data = table2array(atmo(:,5)); % m/s
% Cp_data = table2array(atmo(:,6)); % J/kg/K
% Visc_data = table2array(atmo(:,7)); % Ns/m^2

% sample accoroding to model z array
T = interp1(z_data,T_data,z_array,'spline');
rho = interp1(z_data,rho_data,z_array,'spline');
p = interp1(z_data,P_data,z_array,'spline');
U = interp1(z_data,U_data,z_array,'spline');
% Cp = interp1(z_data,Cp_data,z_array,'spline');
% Visc = interp1(z_data,Visc_data,z_array,'spline');

% rough estimates:
C = 200.*ones(size(z_array));
R = 191.8.*ones(size(z_array));
gamma = 1.29.*ones(size(z_array));
H = 11.1.*ones(size(z_array));
kvisc = 10000.*ones(size(z_array)); % ignore, only used for sponge layer
thermdiffus = kvisc; % ignore
