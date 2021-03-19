function[time,domain,dCFL,difCFL,flags,atmo,wind,forcing] = config()

%% ---- Time settings ----
time.Tmin  = 0;    % Initial time
time.Tmax  = 8000; % Final time in seconds
time.skipT = 30;  % Number of seconds to skip storing results
% computation happens after every dt but only limited data is stored
time.n = 0;       % First Step n = 0 (n counts the current stored frame)
time.t = time.Tmin;    % First time t = Tmin
time.nframe = 1;  % First frame = 1 (nframe is the total no of stored frames)
time.T_arr(time.nframe) = 0; % T_arr is the time array of stored frames

%% ---- Domain settings ----
domain.Xmin = -20000;
domain.Xmax = 20000;
domain.Zmin = 0;
domain.Zmax = 160000;
domain.dx = 500; % horizontal resolution
domain.dz = 500; % vertical resolution
domain.SpongeHeight = 50000; % sponge layer thickness in meters

%% ---- CFL ----
dCFL = 0.8; % desired Courant-Friedrichs-Lewy number
difCFL = 0.2; % CFL for diffusion problem

%% ---- Viscous phenomenon flags ----
flags.IsTopSpongeLayer = 1; % flag to include a 50 km sponge layer on top
flags.IsViscosity = 0;% flag to solve for molecular viscosity
flags.IsConduction = 0; % flag to solve for thermal conduction  

%% ---- Background atmosphere ----
% Using Earth isothermal model
[~,atmo.rho0,atmo.P0,atmo.R,atmo.gamma,atmo.kinvisc,atmo.H,atmo.C] = Earth_isothermal(Z);
 
%% ---- Background wind ----
% Gaussian wind shear
wind.u_max = 0;    % wind amplitude (m/s) 
wind.u_zloc = 100000;    % z location of wind peak (m)
wind.u_sig = 10000;    % stdev of wind profile (m)
wind.wind = u_max.*exp(-(Z-u_zloc).^2./(2*u_sig^2));    % also a matrix of size X=Z

%% ---- Wave forcing ----
% A lower boundary Source is simulated as Gaussian w perturbation
forcing.no = false;     %if true, no forcing is applied
forcing.amp = 0.001;      % amplitude (m/s)
forcing.omega = 0.007;  % centered frequency
kx = 2*pi / (Xmax-Xmin);    % One horizontal wavelength per domain is set (lambda_x = x domain length)
forcing.kxx = x_c.*kx;  % computing kx*x
forcing.t0 = 1200;      % time at forcing maxima (s)
forcing.sigmat=600;     % forcing half width time (s)
