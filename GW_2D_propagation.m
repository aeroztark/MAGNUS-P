
clc
clear all
close all

% import appropriate simulation configuration file
config_baseline;

%% Additional configuration from inputs of the config file

% note on indexing: X(row,col) --> X(z_ind, x_ind)

% managing z axis indices
LastDomainZindex = find(z_c > ZDomainEnd & z_c < ZDomainEnd+dz)-1; % last Z index for the 'physically valid' domain
FirstSpongeZindex = LastDomainZindex + 1;
if IsTopSpongeLayer == 0  % if no sponge layer is implemented, just take last 2 indices out for top BCs
    LastDomainZindex = LastDomainZindex - 2;
end

% setting viscosity coefficient constant in sponge layer to prevent diffusion timestep being too low
% may not be a concern with implicit method
 if IsTopSpongeLayer == 1
     kinvisc(FirstSpongeZindex:end,:) = kinvisc(LastDomainZindex,1);
 end

% model gravity to maintain hydrostatic equilibrium initially (g dimension is Z-1 x X-1)
g = (P0(2:end,1)-P0(1:end-1,1))./(-0.5*dz*(rho0(2:end,1)+rho0(1:end-1,1)));
g = repmat(g,1,size(X,2)-1);

% ---- initial timestep calculation
dt = dCFL.*min(dx,dz)./max(max(C));   %limited by speed of sound

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
    S=Source(0.5*(Q(jD,iD,:)+Q(jD+1,iD,:)),g(jD,iD),X,Z,t); % compute source
    % half step in z
    Qs(jD,iD,:)=0.5*(Q(jD,iD,:)+Q(jD+1,iD,:))-(dt/(2*dz))*(G(jD+1,iD,:)-G(jD,iD,:))+(dt/2)*S(:,:,:);
    G=Gflux(Qs);    % update flux
    S=Source(Qs,g,X,Z,t); % update source
    % full step in z
    Q(jD,iD,:)=Q(jD,iD,:)-(dt/dz)*(G(jD,iD,:)-G(jD-1,iD,:))+dt*0.5*(S(jD,iD,:)+S(jD-1,iD,:));
    % apply BCs
    Q=bc(Q,t);
    
    % Solve for diffusion terms
    % Molecular Diffusion
    if IsViscosity == 1
        Q = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit);
    end
    
    % Thermal conductivity
    if IsConduction == 1
        Q = ThermalConduction(thermdiffus,difCFL,T0,dt,dx,dz,x_c,z_c,jD,iD,Q,t,IsDiffusionImplicit);
    end
    
    % Sponge layer implementation
    if IsTopSpongeLayer == 1
        Q = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit);        
    end
    
    % ---- Update time ----
    t=t+dt;
    n=n+1;
    
    % Update advective (main loop) timestep (adaptive)
    dt = dCFL.*min(dx,dz)./(max(max(C))+max(max(max(abs(Q(:,:,2:3)./Q(:,:,1)))))); % when Q(:,:,1) -> 0, dt becomes 0 & program breaks
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
% all outputs are only taken from indices (3:end-2) in X and (3:LastDomainZindex) in Z since that is the
% computational domain, after excluding 2 ghost cells on either sides and sponge layer, if implemented.

% These values are 3D arrays (z-x-t)
KE = squeeze(0.5*(Q_save(3:LastDomainZindex,3:end-2,2,:).^2+Q_save(3:LastDomainZindex,3:end-2,3,:).^2)./Q_save(3:LastDomainZindex,3:end-2,1,:));
P_PERT = (squeeze(Q_save(3:LastDomainZindex,3:end-2,4,:))-KE).*(gamma(3:LastDomainZindex,3:end-2)-1)-P0(3:LastDomainZindex,3:end-2);
T_PERT = P_PERT./(R(3:LastDomainZindex,3:end-2).*squeeze(Q_save(3:LastDomainZindex,3:end-2,1,:)));
U = squeeze(Q_save(3:LastDomainZindex,3:end-2,2,:)./Q_save(3:LastDomainZindex,3:end-2,1,:));    %horiz. wind (not perturbation)
W = squeeze(Q_save(3:LastDomainZindex,3:end-2,3,:)./Q_save(3:LastDomainZindex,3:end-2,1,:));    % vertical wind (not perturbation)

SCALING_FACTOR = sqrt(rho0(3:LastDomainZindex,3:end-2)./rho0(3,3:end-2)); % an 2d Z-X matrix
Z_KM = z_c(3:LastDomainZindex)./1000; % grid center arrays for plotting the computational domain
X_KM = x_c(3:end-2)./1000;

%% Optional plots
figure % wave travel plot for x in the middle of domain
contourf(T_arr,Z_KM,squeeze(W(:,length(X_KM)/2,:)).*SCALING_FACTOR(:,length(X_KM/2)),50,'Edgecolor','none')
xlabel('time (s)')
ylabel('z (km)')

% BVf = BV(g(:,1),C(:,1),gamma(:,1));
% Ri = RichardsonNumber(BVf,LastDomainZindex,U(:,length(X_KM)/2,500),500);
% figure
% plot(Ri,Z_KM,'LineWidth',2)
% xlim([-1 1])

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
    if ~(forcing.verticalvelocity)   % i.e. if no vertical velocity forcing, use reflective BC for rho*w at domain bottom
        Q(1:2,:,3) = -Q(3,:,3).*(rho0(1:2,:)./rho0(3,:)).^(0.5); 
    else % enforce vertical velocity forcing
        w = forcing.amp.*cos(forcing.omega.*(t-forcing.t0)-forcing.kxx).*exp(-(t-forcing.t0)^2./(2*forcing.sigmat^2));
        %w = Tsunami_forcing(t); 
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
function S = Source(Q,g,X,Z,t)
    global forcing
    S(:,:,1) = 0.*Q(:,:,1);
    S(:,:,2) = 0.*Q(:,:,1);
    S(:,:,3) = -Q(:,:,1).*g;
    if ~(forcing.thermal) %if thermal forcing is not applied 
        S(:,:,4) = -Q(:,:,3).*g;   % no thermal forcing (just -rho*g*w)
    else
        [a,b,~] = size(Q);
        S(:,:,4) = -Q(:,:,3).*g + Q(:,:,1).*forcing.amp.*exp(-((X(1:a,1:b)-forcing.x0).^2)./(2*forcing.sigmax^2)).*exp(-((Z(1:a,1:b)-forcing.z0).^2)./(2*forcing.sigmaz^2)).*exp(-((t-forcing.t0).^2)./(2*forcing.sigmat^2));
    end
end

%% ---- Viscous terms ----

function[Q] = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit)
global wind
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
    max_visc = max(max(kinvisc));  %max value of viscosity in the domain
    
    if IsDiffusionImplicit == 1
        N_substeps = ceil(dt/10); % run 10 substeps for every advection timestep
    else
        N_substeps = ceil(dt*max_visc/(difCFL*min(dx,dz)^2));   %no of substeps required to solve diffusion equation...
        ... based on Von Neumann Number (dt = N_substeps x dt_sub)
    end

    % Main Substepping Loop
    for m = 1:N_substeps
        
        % Substep Timestep
        dt_sub = dt/N_substeps;
        
            % x-split for diffusion equation ----
            % compute intermediate values for ease:
            Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1); % compute P/(gamma-1) 
            Q(:,:,2:3) = Q(:,:,2:3)./Q(:,:,1);  % compute velocity
            u_diff = Q(:,:,2) - wind;
            w_diff = Q(:,:,3);
            if IsDiffusionImplicit == 1
                % Using direct Implicit method:
                % first create coefficient matrix
                A = CreateImplicitMatrix(Q(:,:,2),jD,iD,kinvisc.*(dt_sub/(dx^2)));
                % then use LU decomposition to solve the linear system
                Q(jD,iD,2) = wind(jD,iD) + A \ u_diff(jD,iD);
                Q(jD,iD,3) = A \ w_diff(jD,iD);
            else
                % Using an explicit scheme: forward Euler in time and centrered difference in space
                Q(jD,iD,2:3) = Q(jD,iD,1).*(Q(jD,iD,2:3)+kinvisc(jD,iD).*(dt_sub/dx^2).*(Q(jD,iD+1,2:3)-2.*Q(jD,iD,2:3)+Q(jD,iD-1,2:3))); % get updated rho*u and rho*w
            end
            
             Q(jD,iD,4) = Q(jD,iD,4)+0.5*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1); % get updated E
            
            % apply BCs
            Q = bc(Q,t);

            % z-split for diffusion equation ----
            Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1);
            Q(:,:,2:3)=Q(:,:,2:3)./Q(:,:,1);
            u_diff = Q(:,:,2) - wind;
            w_diff = Q(:,:,3);
             if IsDiffusionImplicit == 1
                % Using direct Implicit method:
                % first create coefficient matrix
                A = CreateImplicitMatrix(Q(:,:,2),jD,iD,kinvisc.*(dt_sub/(dz^2)));
                % then use LU decomposition to solve the linear system
                Q(jD,iD,2) = wind(jD,iD) + A \ u_diff(jD,iD);
                Q(jD,iD,3) = A \ w_diff(jD,iD);
             else
                 % Explicit method
                Q(jD,iD,2:3) = Q(jD,iD,1).*(Q(jD,iD,2:3)+kinvisc(jD,iD).*(dt_sub/dz^2).*(Q(jD+1,iD,2:3)-2.*Q(jD,iD,2:3)+Q(jD-1,iD,2:3)));
             end
             
             Q(jD,iD,4) = Q(jD,iD,4)+0.5*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1);
         
        % apply BCs
        Q = bc(Q,t);
    
    end
end

function[Q] = ThermalConduction(thermdiffus,difCFL,T_ref,dt,dx,dz,x_c,z_c,jD,iD,Q,t,IsDiffusionImplicit)

global R gamma P0
    % This function solves the diffusion equation for thermal conduction
    % inputs: thermdiffus -> array of thermal diffusivity values
             % T_ref -> backgrounf temperature field
             % difCFL -> CFL number for solving diffusive equations (must be <0.5)
             % dt     -> time step for the main (advective) loop
             % dz     -> spatial grid resolution
             % jD     -> indices for computation in the grid
             % Q      -> Array of prognostic vars at a timestep
             % t      -> current simulation epoch time
     % outputs: Q-> updated Q after solving for viscosity and updating u and E   
     
    % First calculating the number of sub-steps for integration of diffusion equation 
    max_diffusivity = max(max(thermdiffus));  %max value of viscosity in the domain
    
    if IsDiffusionImplicit == 1
        N_substeps = ceil(dt/10); % run 10 substeps for every advection timestep
    else
        N_substeps = ceil(dt*max_diffusivity/(difCFL*min(dx,dz)^2));   %no of substeps required to solve diffusion equation...
        ... based on Von Neumann Number (dt = N_substeps x dt_sub)
    end

    % Main Substepping Loop
    for m = 1:N_substeps
        
        % Substep Timestep
        dt_sub = dt/N_substeps;
        
            % x-split for diffusion equation ----
            % compute intermediate values for ease:
            Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1); % compute P/(gamma-1) 
            T = (Q(:,:,4).*(gamma-1))./(R.*Q(:,:,1));  % compute T
            T_diff = T - T_ref; % diffusion will be applied only to deviation from reference state
            
            if IsDiffusionImplicit == 1
                % Using direct Implicit method:
                A = CreateImplicitMatrix(T_diff,jD,iD,thermdiffus.*(dt_sub/(dx^2)));
                % then use LU decomposition to solve the linear system
                T_diff(jD,iD) = A \ T_diff(jD,iD);
            else
                % Using an explicit scheme: forward Euler in time and centrered difference in space
                T_diff(jD,iD) = T_diff(jD,iD) + thermdiffus(jD,iD).*(dt_sub/dx^2).*(T_diff(jD,iD+1)-2.*T_diff(jD,iD)+T_diff(jD,iD-1)); % get updated T
            end
            
            P = Q(:,:,1).*R.*(T_diff + T_ref); % get updated P
            Q(jD,iD,4) = P(jD,iD)./(gamma(jD,iD)-1) + 0.5.*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1); % get updated E

            % apply BCs
            Q = bc(Q,t);

            % z-split for diffusion equation ----
            Q(:,:,4) = Q(:,:,4)-0.5*(Q(:,:,2).^2+Q(:,:,3).^2)./Q(:,:,1); % compute P/(gamma-1) 
            T = (Q(:,:,4).*(gamma-1))./(R.*Q(:,:,1));  % compute T
            T_diff = T - T_ref; % diffusion will be applied only to deviation from reference state
            
            if IsDiffusionImplicit == 1
                % Using direct Implicit method:
                A = CreateImplicitMatrix(T_diff,jD,iD,thermdiffus.*(dt_sub/(dz^2)));
                % then use LU decomposition to solve the linear system
                T_diff(jD,iD) = A \ T_diff(jD,iD);
            else
                % Using an explicit scheme: forward Euler in time and centrered difference in space
                T_diff(jD,iD) = T_diff(jD,iD) + thermdiffus(jD,iD).*(dt_sub/dz^2).*(T_diff(jD+1,iD)-2.*T_diff(jD,iD)+T_diff(jD-1,iD)); % get updated T
            end
            P = Q(:,:,1).*R.*(T_diff + T_ref); % get updated P
            Q(jD,iD,4) = P(jD,iD)./(gamma(jD,iD)-1) + 0.5*(Q(jD,iD,2).^2+Q(jD,iD,3).^2)./Q(jD,iD,1); % get updated E
        end
   
        % apply BCs
        Q = bc(Q,t);

    end


function[A] = CreateImplicitMatrix(Q,jD,iD,F)
% function to create the coefficient matrix for direct implicit solution of
% 1D diffusion equation (Dimensionally split diffusion equation will be solved using LU factorization)
% Inputs: Q -> Vector of quantities
%         jD,iD -> domain indices
%         F - mesh Fourier number = K*dt/(ds^2)

[nrows,ncols] = size(Q(jD,iD));
% Set the main (zero-th) diagonal as 2F+1
A = spdiags(1+2.*F(:,1),0,zeros(nrows));
% Set the first (+1th and -1th) diagonals as -F
A = spdiags(-F(:,1),1,A);
A = spdiags(-F(:,1),-1,A);
end

