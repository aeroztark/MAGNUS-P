% Config file for Mars cloud case
%% ---- Time settings ----
Tmin  = 0;    % Initial time
Tmax  = 18000; % Final time in seconds (3 hrs)
skipT = 180;  % Number of seconds to skip storing results (1 min)
% computation happens after every dt but only limited data is stored
n = 0;       % First Step n = 0 (n counts the current stored frame)
t = Tmin;    % First time t = Tmin
nframe = 1;  % First frame = 1 (nframe is the total no of stored frames)
T_arr(nframe) = 0; % T_arr is the time array of stored frames

%% ---- Viscous phenomenon flags ----
IsTopSpongeLayer = 1; % flag to include a sponge layer on top
IsViscosity = 0;% flag to solve for molecular viscosity
IsConduction = 0; % flag to solve for thermal conduction  
IsDiffusionImplicit = 0;
%% ---- Domain settings ----
% note on indexing: X(row,col) --> X(z_ind, x_ind)
Xmin = 0;
Xmax = 550000;
Zmin = 0;
Zmax = 15000;
dx = 500; % horizontal resolution
dz = 250; % vertical resolution
SpongeHeight = 20000; % sponge layer thickness in meters

ZDomainEnd = Zmax; % last value of Z for ' physically valid' domain (pre-sponge)

if IsTopSpongeLayer == 1
    Zmax = Zmax + SpongeHeight;    % extend Z by 50 km for computations
end

Xdomain = Xmax-Xmin;
Zdomain = Zmax-Zmin;
% Define positions at cell centers
x_c = Xmin-3*dx/2:dx:Xmax+3*dx/2;   %2 centre points outside on either boundary
z_c = Zmin-3*dz/2:dz:Zmax+3*dz/2;
[X,Z] = meshgrid(x_c,z_c);    % grid of cell centers
[J,I] = size(X);  %J gives no of z levels and I gives no of x levels

%% ---- CFL ----
dCFL = 0.85; % desired Courant-Friedrichs-Lewy number
difCFL = 0.4; % CFL for diffusion problem

%% ---- Background atmosphere ----
global g R P0 rho0 gamma C; 

% If using Earth isothermal model
% [T0,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C] = Earth_isothermal(Z);
 
% Using Mars model from MCD
  [~,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C,U0] = Mars_MOLApass260(Z);
%  [~,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C,U0] = Mars_MOLApass390(Z);
  
 T0 = 1.01.*(90 + ((104.5-Z./1000)./1.853)); % 2% subsaturated
 
%% ---- Background wind ----
% only horizontal wind is specified -> time invariant. Vertical wind is zero.
% will need recheck equations to subtract from rho*u and rho*w if this changes 

global wind
% Gaussian wind shear
% u_max = 0;    % wind amplitude (m/s) 
% u_zloc = 100000;    % z location of wind peak (m)
% u_sig = 10000;    % stdev of wind profile (m)
% wind = u_max.*exp(-(Z-u_zloc).^2./(2*u_sig^2));    % wind profile, also a matrix of size X=Z

% linear wind shear
%wind = -U0;
wind = 10.*ones(size(U0));

%% ---- Wave forcing ----
% A tsunami forcing function is called in the main file
global forcing
forcing.no = false;     %if true -> no forcing is applied

% idealized ridge trough
% H = 0.5; % km
% L = 7; %km 
% hq = -H.*exp(-((x_c./1000).^2)./(L^2));

%load('Smooth_pass260_topo.mat'); % load smoothened topography
% sample topography at model x scale
%hq = interp1(x,h,x_c./1000,'pchip'); % for 260

% load('Smooth_pass390_topo.mat');
% h1 = interp1(x,h,0:3:200,'pchip'); % for 390 (2-stage snoothening to remove hiccupy gradient)
% hq = interp1(0:3:200,h1,x_c./1000,'pchip'); % for 390
%hq = interp1(x,h,x_c./1000,'pchip');

load('topo207.mat'); % load smoothened topography
% sample topography at model x scale
hq = interp1(x,h,x_c./1000,'pchip'); % for 260

% Forcing from the obtained topography
dh_dx = diff(hq)./(dx/1000); % gradient
dh_dx = [dh_dx(1), dh_dx]; % append the same at position 1 to make the same length vector

forcing.topo_w = 10.*dh_dx;  % surface wind x topography horizontal gradient
