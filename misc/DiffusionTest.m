% This script compares the different numerical methods against the exact
% solution of a diffusion problem

clear

% Diffusion testing
k = 100;
Lx = 100;
Ly = 100;
dx = 1;
dy = 1;
param = -(k*pi^2)/(Lx^2 + Ly^2); % for ease of computation
x = 1:dx:Lx;
y = 1:dy:Ly;
t = 0;
dt = 1;
T_end = 100;

t_array = 0:dt:T_end;
max_err = zeros(size(t_array));
rms_err = zeros(size(t_array));

u_SOR = zeros(length(x),length(y));
u_exact = zeros(length(x),length(y));
u_explicit = zeros(length(x),length(y));

rho_para = (cos(pi/length(x)) + ((dx/dy)^2)*cos(pi/length(y)))/(1+(dx/dy)^2);
omega_optimal = 2/(1+sqrt(1-rho_para^2));

% IC
for i = 1:Lx
    for j = 1:Ly
        u_exact(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        u_SOR(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        u_explicit(i,j) = 10.*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
    end
end

%% time loop
for ind = 1:length(t_array)
    
    % BCs
    u_exact(1,:) = 0;
    u_exact(end,:) = 0;
    u_exact(:,1) = 0;
    u_exact(:,end) = 0;
      
%     % exact solution
    for i = 1:Lx
        for j = 1:Ly
            u_exact(i,j) = 10.*exp(param*t).*sin((pi/Lx).*x(i)).*sin((pi/Ly).*y(j)) ;
        end
    end

%u = Implicit_Diffusion(length(x),length(y),dx,dy,dt,k,u);
  
% SOR method:
    u_SOR(1,:) = 0;
    u_SOR(end,:) = 0;
    u_SOR(:,1) = 0;
    u_SOR(:,end) = 0;
    
    % SOR solution
    %u_SOR = SOR(u_SOR,dt,dx,dy,x,y,k,1e-6,1);
    u_SOR = gauss_seidel(dx,dy,dt,2:length(x)-1,2:length(y)-1,u_SOR,k,1e-6,1); % gauss seidel

    SORdiff = (u_exact - u_SOR);
    SORmax_err(ind) = max(max(abs(SORdiff)));
    SORrms_err(ind) = rms(reshape(SORdiff,[length(x)*length(y),1]));

% Explicit method
  u_explicit = ExplicitMethod(k,0.5,u_explicit,dt,dx,dy,2:length(y)-1,2:length(x)-1);
    Explicitdiff = (u_exact - u_explicit);
    Explicitmax_err(ind) = max(max(abs(Explicitdiff)));
    Explicitrms_err(ind) = rms(reshape(Explicitdiff,[length(x)*length(y),1]));
%     surf(u_SOR)
%     zlim([0 10])
% %     contourf(u,50,'Edgecolor','none')
% %     caxis([0 1])
%     colorbar
%      pause(0.01)
ind = ind+1;
t = t+dt;
end


figure
subplot(1,2,1)
plot(t_array,SORmax_err)
hold on
plot(t_array,Explicitmax_err)
legend('SOR','Explicit')
subplot(1,2,2)
plot(t_array,SORrms_err)
hold on
plot(t_array,Explicitrms_err)
legend('SOR','Explicit')
%% ----------------------------
function[u_prev_t] = SOR(u_old,dt,dx,dy,x,y,alpha,tolerance,omega)

% This function does not work as expected for some reason
% inputs: u_old -> variable matrix

Fx = alpha.*(dt/dx^2); % mesh Fourier numbers
Fy = alpha.*(dt/dy^2);
k1 = alpha*dt/(dx^2);
k2 = alpha*dt/(dy^2);

error = 1; % init

u_prev_t = u_old;
u_prev_iter = u_old;
u_star = zeros(size(u_old));
u_new_iter = zeros(size(u_old));

while error > tolerance
    
    for i = 2:(length(x)-1) % only evaluate at interior points
        for j = 2:(length(y)-1)
            u_star(i,j) = (u_prev_t(i,j) + Fx.*(u_old(i-1,j) + u_prev_iter(i+1,j))+...
                        Fy.*(u_prev_t(i,j-1) + u_prev_iter(i,j+1)))./(1+2.*Fx + 2.*Fy);              
             u_new_iter = omega*u_star + (1-omega)*u_prev_t;
        end
    end
    
    % update value for next iteration (same time)
    u_prev_iter = u_new_iter;
    error = max(max(abs(u_new_iter - u_prev_iter)));
end
% when error < tolerance, output the latest iteration value as the result
% of this time step
u_prev_t = u_new_iter;
end

%% ----------------------
function[u] = ExplicitMethod(k,difCFL,u,dt,dx,dy,jD,iD)


    % First calculating the number of sub-steps for integration of diffusion equation 
    N_substeps = ceil(dt*k/(difCFL*min(dx,dy)^2));   %no of substeps required to solve diffusion equation...
        ... based on Von Neumann Number (dt = N_substeps x dt_sub)
    
    % Main Substepping Loop
    for m = 1:N_substeps
        
        % Substep Timestep
        dt_sub = dt/N_substeps;
        
        % x-split for diffusion equation ----

        % Using an explicit scheme: forward Euler in time and centrered difference in space
        u(jD,iD) = u(jD,iD) + k.*(dt_sub/dx^2).*(u(jD,iD+1)-2.*u(jD,iD)+u(jD,iD-1)); % get updated T
       
        % apply BCs
        u(1,:) = 0;
        u(end,:) = 0;
        u(:,1) = 0;
        u(:,end) = 0;
        
        % z-split for diffusion equation ----
        u(jD,iD) = u(jD,iD) + k.*(dt_sub/dy^2).*(u(jD+1,iD)-2.*u(jD,iD)+u(jD-1,iD)); % get updated T
          

        % apply BCs
        u(1,:) = 0;
        u(end,:) = 0;
        u(:,1) = 0;
        u(:,end) = 0;
    
    end
end