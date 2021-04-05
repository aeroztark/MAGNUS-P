    
function[U_diff] = gauss_seidel(dx,dz,dt,iD,jD,U_diff,alpha,tolerance,omega)

% Function to solve 2D diffusion equation implicitly using
% Gauss-Seidel iterative solver: to be used for solution of molecular
% viscosity and thermal conduction equations
%   dU/dt = alpha*[Uxx + Uzz

% inputs:
% dx,dz,dt -> space and time steps
% iD,jD -> domain indices from the simulation script
% U_diff -> quantity to be solved for (prognostic var)
% alpha -> viscosity coefficient (kinematic visc or themal diffusivity)
% tolerance -> suitable tolerance for the solution 
% omega  -> relaxation factor (set as 1 for Gauss Seidel)

     
    U_old = U_diff;                  % for updation old values in convergence loop
    U_prev = U_diff;
    n_iter = 1;            % to count the number of total iterations            

    % Mesh Fourier numbers
    Fx = alpha.*(dt/dx^2);         % for ease of calculation
    Fz = alpha.*(dt/dz^2);
    
   
    error = 1;              % error initialized to 1 before each time loop
    
    % convergence loop
    while error > tolerance
        
    % nodal loop
    % implicit scheme using Gauss-Seidel
        U_diff(jD,iD) = (U_prev(jD,iD) + Fx(jD,iD).*(U_diff(jD,iD-1) + U_old(jD,iD+1))+...
                        Fz(jD,iD).*(U_diff(jD-1,iD) + U_old(jD+1,iD)))./(1+2.*Fx(jD,iD) + 2.*Fz(jD,iD));
        % if doing SOR:         
        %U_diff(jD,iD) = omega*U_diff(jD,iD) + (1-omega)*U_prev(jD,iD);    
        
        % convergence criterion
        error = max(max(abs(U_diff - U_old)));
        
        % updating old values
        U_old = U_diff;
        n_iter = n_iter +1;
    end
 
end