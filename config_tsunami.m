% Config file for Tsunami case

%% ---- Time settings ----
Tmin  = 0;    % Initial time
Tmax  = 10800; % Final time in seconds (3 hrs)
skipT = 60;  % Number of seconds to skip storing results (1 min)
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
Xmin = -3000000;
Xmax = 3000000;
Zmin = 0;
Zmax = 200000;
dx = 10000; % horizontal resolution
dz = 500; % vertical resolution
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
dCFL = 0.8; % desired Courant-Friedrichs-Lewy number
difCFL = 0.2; % CFL for diffusion problem

%% ---- Background atmosphere ----
global g R P0 rho0 gamma C; 

% If using Earth isothermal model
% [T0,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C] = Earth_isothermal(Z);
 
% If using Earth MSIS model
 [T0,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C,~,~] = Earth_MSIS(Z,3.5,96,2004,361,7200);
 
%% ---- Background wind ----
% only horizontal wind is specified -> time invariant. Vertical wind is zero.
% will need recheck equations to subtract from rho*u and rho*w if this changes 

global wind
% Gaussian wind shear
u_max = 0;    % wind amplitude (m/s) 
u_zloc = 100000;    % z location of wind peak (m)
u_sig = 10000;    % stdev of wind profile (m)
wind = u_max.*exp(-(Z-u_zloc).^2./(2*u_sig^2));    % wind profile, also a matrix of size X=Z

% linear wind shear
% wind = linspace(0,u_max,length(z_c));
% wind = repmat(wind',1,length(x_c));    % also a matrix of size X=Z

%% ---- Wave forcing ----
% A tsunami forcing function is called in the main file
global forcing
forcing.no = false;     %if true -> no forcing is applied
% forcing.amp = 0.001;      % amplitude (m/s)
% forcing.omega = 0.007;  % centered frequency
% kx = 2*pi / (Xmax-Xmin);    % One horizontal wavelength per domain is set (lambda_x = x domain length)
% forcing.kxx = x_c.*kx;  % computing kx*x
% forcing.t0 = 1200;      % time at forcing maxima (s)
% forcing.sigmat=600;     % forcing half width time (s)
