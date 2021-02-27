
clc
clear
close all

%% Simulation inputs
% ---- Domain size (in meters) ----
Xmin = -20000;
Xmax = 20000;
Zmin = 0;
Zmax = 160000;
dx = 500; % horizontal resolution
dz = 500; % vertical resolution
% note on indexing: X(row,col) --> X(z_ind, x_ind)

Xdomain = Xmax-Xmin;
Zdomain = Zmax-Zmin;
% Define positions at cell centers
x_c = Xmin-3*dx/2:dx:Xmax+3*dx/2;   %2 centre points outside on either boundary
z_c = Zmin-3*dz/2:dz:Zmax+3*dz/2;
[X,Z] = meshgrid(x_c,z_c);    % grid of cell centers
[J,I] = size(X);  %J gives no of z levels and I gives no of x levels

% ---- Time parameters ----
Tmin  = 0;    % Initial time
Tmax  = 3000; % Final time in seconds
skipT = 30;  % Number of seconds to skip storing results
% computation happens after every dt but only limited data is stored
n = 0;       % First Step n = 0 (n counts the current stored frame)
t = Tmin;    % First time t = Tmin
nframe = 1;  % First frame = 1 (nframe is the total no of stored frames)
T_arr(nframe)=0; % T_arr is the time array of stored frames

% ---- Background atmosphere

global g R P0 rho0 gamma C; 

% Using Earth isothermal model
 [~,rho0,P0,R,gamma,kinvisc,H,C] = Earth_isothermal(Z);
% Using Earth MSIS model
%[~,rho0,P0,R,gamma,kinvisc,H,C,] = Earth_MSIS(Z,10,180,2020,1,0);

% model gravity to maintain hydrostatic equilibrium initially (g dimension is Z-1 x X-1)
g = (P0(2:end,1)-P0(1:end-1,1))./(-0.5*dz*(rho0(2:end,1)+rho0(1:end-1,1)));
g = repmat(g,1,size(X,2)-1);

% Viscosity and Thermal Diffusivity Profiles
% dynvisc = 1.3e-5; % sample molecular viscosity
% Pr = 0.7; % Prandtl number
% kinvisc = dynvisc./rho0;
% tdiffus = kinvisc./Pr;

IsViscosity = 0;% flag to solve for molecular viscosity
IsConduction = 0; % flag to solve for thermal conduction   

% ---- Background wind (modelled as horizontal wind with vertical Gaussian or Linear profile)
global wind

% Gaussian wind shear
 u_max = 100;    % wind amplitude (m/s) 
u_zloc = 100000;    % z location of wind peak (m)
u_sig = 10000;    % stdev of wind profile (m)
wind= u_max.*exp(-(Z-u_zloc).^2./(2*u_sig^2));    % also a matrix of size X=Z

% linear wind shear
% wind = linspace(0,u_max,length(z_c));
% wind = repmat(wind',1,length(x_c));

% ---- Wave forcing
global forcing
% A lower boundary Source is simulated as Gaussian w perturbation
forcing.no = false;     %if true, no forcing is applied
forcing.amp = 0.002;      % amplitude (m/s)
forcing.omega = 0.012;  % centered frequency
kx = 2*pi / (Xmax-Xmin);    % One horizontal wavelength per domain is set (lambda_x = x domain length)
forcing.kxx = x_c.*kx;  % computing kx*x
forcing.t0 = 1200;      % time at forcing maxima (s)
forcing.sigmat=600;     % forcing half width time (s)

% ---- timestep and CFL ----
dCFL = 0.8; % desired Courant-Friedrichs-Lewy number
difCFL = 0.2; % CFL for diffusion problem
dt = dCFL.*min(dx,dz)./max(C,[],'all');   %limited by speed of sound

%------End of Inputs--------------------------- 

%% Simulation setup

% Initializing arrays:
    %our PDE system is: dQ/dt + dF/dx + dG/dz = S 
    % all these are zero arrays of size X in 4D 
F = cat(3,0.*X,0.*X,0.*X,0.*X);   %Fluxes for x-split(concantenate along 3rd dim to yield 4D array)
G = F;                            %Fluxes for y-split
Q = F;                            %Solution Variables (rho, rho*u, rho*w, E)
S = F;                            %Source Terms
Q_save = zeros(J,I,4,nframe);     %4D array to hold Q results at certain cadence

% ---- Initial Conditions ----
P_pert=0.*X;    % zero pressure perturbation
%P_pert=1000*exp(-.5*((X-0).^2)./(Xdomain/10).^2).*exp(-.5*((Z-10000./2).^2)./(Zdomain/10).^2);
Q(:,:,1) = rho0;    %rho
Q(:,:,2) = rho0.*wind;    % rho*u
Q(:,:,3) = 0;             % rho*w (forcing is added in BCs)
Q(:,:,4) =(P_pert+P0.*X.^0)./(gamma-1)+0.5*rho0.*wind.^2; % E for ideal gas

% Set BCs for first time
Q = bc(Q,0);

% Set domain indices for evaluations (leave out 1st and last gridcenter)
iD = 2:I-1;   % x indices
jD = 2:J-1;   % z indices

% Store initial state
Q_save(:,:,:,nframe) = Q;   %currently nframe is 1

%% Computations (Lax-Wendroff 2-step)

while t < Tmax
    
    % ---- x-split ---- (no source used in x split since our sources are height dependent)
    F=Fflux(Q); % compute flux F    
    % half-step
    Qs(jD,iD,:)=0.5*(Q(jD,iD,:)+Q(jD,iD+1,:))-(dt/(2*dx))*(F(jD,iD+1,:)-F(jD,iD,:));
    F=Fflux(Qs); % update flux
    % full-step in x
    Q(jD,iD,:)=Q(jD,iD,:)-(dt/dx)*(F(jD,iD,:)-F(jD,iD-1,:));
    % apply BCs
    Q=bc(Q,t);
    
    % z-split
    G=Gflux(Q); % compute flux G
    S=Source(0.5*(Q(jD,iD,:)+Q(jD+1,iD,:)),g(jD,iD)); % compute source
    % half step in z
    Qs(jD,iD,:)=0.5*(Q(jD,iD,:)+Q(jD+1,iD,:))-(dt/(2*dz))*(G(jD+1,iD,:)-G(jD,iD,:))+(dt/2)*S(:,:,:);
    G=Gflux(Qs);    % update flux
    S=Source(Qs,g); % update source
    % full step in z
    Q(jD,iD,:)=Q(jD,iD,:)-(dt/dz)*(G(jD,iD,:)-G(jD-1,iD,:))+dt*0.5*(S(jD,iD,:)+S(jD-1,iD,:));
    % apply BCs
    Q=bc(Q,t);
    
    % Solve for diffusion terms
    % Molecular Diffusion
    if IsViscosity ~= 0
        Q = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t);
    end
    
    % Thermal conductivity
    
    
    % ---- Update time ----
    t=t+dt;
    n=n+1;
    
    % Update advective (main loop) timestep (adaptive)
    dt = dCFL.*min(dx,dz)./(max(C,[],'all')+max(abs(Q(:,:,2:3)./Q(:,:,1)),[],'all')); % Not Fool-Proof (!)
    if ((isnan(dt)) || (dt==0))
        break;
    end
    
    % Store results
    if (mod(int32(t),int32(skipT))==0) && (int32(T_arr(nframe))~=int32(t))
        nframe=nframe+1;
        Q_save(:,:,:,nframe)=Q;
        T_arr(nframe)=t;   
        disp(['dt=',num2str(dt),'(s); Time Step n=',num2str(n),'; Time t=',num2str(t),'(s)']);
    end
    
end

%% Outputs (sim outputs are in all caps)
% compute fluid properties from saved Q data (rho, rho*u, rho*w and E)
% all outputs are only taken from indices (3:end-2) since that is the
% computational domain, after excluding 2 ghost cells on either sides.

% These values are 3D arrays (z-x-t)
KE = squeeze(0.5*(Q_save(3:end-2,3:end-2,2,:).^2+Q_save(3:end-2,3:end-2,3,:).^2)./Q_save(3:end-2,3:end-2,1,:));
P_PERT = (squeeze(Q_save(3:end-2,3:end-2,4,:))-KE).*(gamma(3:end-2,3:end-2)-1)-P0(3:end-2,3:end-2);
T_PERT = P_PERT./(R(3:end-2,3:end-2).*squeeze(Q_save(3:end-2,3:end-2,1,:)));
U = squeeze(Q_save(3:end-2,3:end-2,2,:)./Q_save(3:end-2,3:end-2,1,:));
W = squeeze(Q_save(3:end-2,3:end-2,3,:)./Q_save(3:end-2,3:end-2,1,:));

SCALING_FACTOR = sqrt(rho0(3:end-2,3:end-2)./rho0(3,3:end-2)); % an 2d Z-X matrix
Z_KM = z_c(3:end-2)./1000; % grid center arrays for plotting the computational domain
X_KM = x_c(3:end-2)./1000;

%% Plotting
figure(1);
%cmocean balance
Nmax=size(T_arr,2);


for n=1:Nmax
    
    %plot results: Pressure
    subplot(2,2,1);
    contourf(X_KM,Z_KM,P_PERT(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['Pressure Perturbation (Pa) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Temperature
    subplot(2,2,2);
    contourf(X_KM,Z_KM,T_PERT(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['Temperature Perturbation (K) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Horizontal Velocity
    subplot(2,2,3);
    contourf(X_KM,Z_KM,U(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['u (m/s) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Vertical Velocity
    subplot(2,2,4);
    contourf(X_KM,Z_KM,W(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['w (m/s) at time t=',num2str(T_arr(n)),'s']);
    
    pause(0.01)
end

%% Function definitions

%% ---- Boundary Conditions ----
function Q = bc(Q,t)
% This function applies boundary conditions to input Q at time t and
% returns Q after modifying it

% refer to Leveque's textbook for details on BCs

    global g R P0 rho0 gamma C; 
    global wind forcing;
    
    % ---- Bottom ----
    % Outflow condition: Zero order extrapolation + scaling to account for
    % stratification (refer Leveque)
    Q(1:2,:,1) = rho0(1:2,:)+(Q(3,:,1)-rho0(3,:)).*(rho0(1:2,:)./rho0(3,:)).^(0.5); % rho
    Q(1:2,:,2) = rho0(1:2,:).*wind(1:2,:)+(Q(3,:,2)-rho0(3,:).*wind(3,:)).*(rho0(1:2,:)./rho0(3,:)).^(0.5); %rho*u
    
    % Adding bottom forcing for rho*w
    if forcing.no   % i.e. if no forcing, use reflective BC for rho*w at domain bottom
        Q(1:2,:,3) = -Q(3,:,3).*(rho0(1:2,:)./rho0(3,:)).^(0.5); 
    else % enforce forcing
        w = forcing.amp.*cos(forcing.omega.*(t-forcing.t0)-forcing.kxx).*exp(-(t-forcing.t0)^2./(2*forcing.sigmat^2));
        Q(1:2,:,3) = w.*rho0(1:2,:);
    end
    % bottom for E
    Q(1:2,:,4) = P0(1:2,:)./(gamma(1:2,:)-1)+(Q(3,:,4)-P0(3,:)./(gamma(3,:)-1)-0.5*rho0(3,:).*wind(3,:).^2).*(rho0(1:2,:)./rho0(3,:)).^(0.5)+0.5*rho0(1:2,:).*wind(1:2,:).^2;
    
    % ---- Top ----
    % open top (outflow) for all 4 quantities
    Q(end-1:end,:,1) = rho0(end-1:end,:)+(Q(end-2,:,1)-rho0(end-2,:)).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
    Q(end-1:end,:,2) = rho0(end-1:end,:).*wind(end-1:end,:)+(Q(end-2,:,2)-rho0(end-2,:).*wind(end-2,:)).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
    Q(end-1:end,:,3) = Q(end-2,:,3).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5);
    Q(end-1:end,:,4) = P0(end-1:end,:)./(gamma(end-1:end,:)-1)+(Q(end-2,:,4)-P0(end-2,:)./(gamma(end-2,:)-1)-0.5*rho0(end-2,:).*wind(end-2,:).^2).*(rho0(end-1:end,:)./rho0(end-2,:)).^(0.5)+0.5*rho0(end-1:end,:).*wind(end-1:end,:).^2;
    
    % ---- Sides ----
    % periodic for both sides (+x and -x)
    % note that indices 1,2,end-1,end are outside computational domain
    Q(:,2,:) = Q(:,end-2,:);
    Q(:,1,:) = Q(:,end-3,:);
    Q(:,end,:)  = Q(:,4,:);
    Q(:,end-1,:)= Q(:,3,:); 

end

%% ---- Flux terms ----
function F = Fflux(Q)
    global g R P0 rho0 gamma C;  
    
    % ensuring that gamma is of same size as Q (when half step Qs is being passed)
    [a,b,~] = size(Q);
    gamm = gamma(1:a,1:b); % gamm is just gamma of the correct size
    
    % compute pressure from ideal gas equation
    P = (gamm-1).*(Q(:,:,4)-(0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1)));
    
    F(:,:,1) = Q(:,:,2);
    F(:,:,2) = Q(:,:,2).*Q(:,:,2)./Q(:,:,1) + P;
    F(:,:,3) = Q(:,:,2).*Q(:,:,3)./Q(:,:,1);
    F(:,:,4) = (Q(:,:,4)+P).*Q(:,:,2)./Q(:,:,1);
end

function G = Gflux(Q)
    global g R P0 rho0 gamma C; 
    
    % ensuring that gamma is of same size as Q (when half step Qs is being passed)
    [a,b,~] = size(Q);
    gamm = gamma(1:a,1:b); % gamm is just gamma of the correct size
    
    P = (gamm-1).*(Q(:,:,4)-(0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1)));
    
    G(:,:,1) = Q(:,:,3);
    G(:,:,2) = Q(:,:,2).*Q(:,:,3)./Q(:,:,1);
    G(:,:,3) = Q(:,:,3).*Q(:,:,3)./Q(:,:,1)+P;
    G(:,:,4) = (Q(:,:,4)+P).*Q(:,:,3)./Q(:,:,1);
end

%% ---- Source term ----
function S = Source(Q,g)
    S(:,:,1) = 0.*Q(:,:,1);
    S(:,:,2) = 0.*Q(:,:,1);
    S(:,:,3) = -Q(:,:,1).*g;
    S(:,:,4) = -Q(:,:,3).*g;
end

%% ---- Viscous terms ----

function[Q] = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t)
    % This function solves the diffusion equation for molecular viscosity
    % inputs: kinvic -> array of kinematic viscosity values
             % difCFL -> CFL number for solving diffusive equations (must be <0.5)
             % dt     -> time step for the main (advective) loop
             % dz     -> spatial grid resolution
             % jD     -> indices for computation in the grid
             % Q      -> Array of prognostic vars at a timestep
             % t      -> current simulation epoch time
     % outputs: Q-> updated Q after solving for viscosity and updating u and E   
     
    % First calculating the number of sub-steps for integration of diffusion equation 
    max_visc = max(kinvisc,[],'all');  %max value of viscosity in the domain
    N_substeps = ceil(dt*max_visc/(difCFL*min(dx,dz)^2));   %no of substeps required to solve diffusion equation...
        ... based on Von Neumann Number (dt = N_substeps x dt_sub)
    
    % Main Substepping Loop
    for m = 1:N_substeps
        
        % Substep Timestep
        dt_sub = dt/N_substeps;
        
        % x-split for diffusion equation ----
        % compute intermediate values for ease:
        Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1); % compute P/(gamma-1) 
        Q(:,:,2:3) = Q(:,:,2:3)./Q(:,:,1);  % compute velocity
        
        % Using an explicit scheme: forward Euler in time and centrered difference in space
        Q(jD,iD,2:3) = Q(jD,iD,1).*(Q(jD,iD,2:3)+kinvisc(jD,iD).*(dt_sub/dx^2).*(Q(jD,iD+1,2:3)-2.*Q(jD,iD,2:3)+Q(jD,iD-1,2:3))); % get updated rho*u and rho*w
        Q(jD,iD,4) = Q(jD,iD,4)+0.5*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1); % get updated E
        
        % apply BCs
        Q = bc(Q,t);
        
        % z-split for diffusion equation ----

        Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1);
        Q(:,:,2:3)=Q(:,:,2:3)./Q(:,:,1);
        
        Q(jD,iD,2:3) = Q(jD,iD,1).*(Q(jD,iD,2:3)+kinvisc(jD,iD).*(dt_sub/dz^2).*(Q(jD+1,iD,2:3)-2.*Q(jD,iD,2:3)+Q(jD-1,iD,2:3)));
        Q(jD,iD,4) = Q(jD,iD,4)+0.5*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1);
        
        % apply BCs
        Q = bc(Q,t);
    
    end
end
