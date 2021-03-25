    
function[T_diff] = gauss_seidel(dx,dz,dt,iD,jD,T_diff,thermdiffus,tolerance,omega)

% Function to solve 2D Transient state heat conduction implicitly using
% Gauss-Seidel iterative solver


     
    T_old = T_diff;                  % for updation old values in convergence loop
    T_prev = T_diff;
    n_iter = 1;            % to count the number of total iterations            

    % Mesh Fourier numbers
    Fx = thermdiffus.*(dt/dx^2);         % for ease of calculation
    Fz = thermdiffus.*(dt/dz^2);
    
   
    error = 1e6;              % error initialized to 1 before each time loop
    
    % convergence loop
    while error > tolerance
    % nodal loop
    % implicit scheme using Gauss-Seidel
        T_diff(jD,iD) = (T_prev(jD,iD) + Fx.*(T_diff(jD,iD-1) + T_old(jD,iD+1))+...
                        Fz.*(T_diff(jD-1,iD) + T_old(jD+1,iD)))./(1+2.*Fx + 2.*Fz);
         T_diff(jD,iD) = omega*T_diff(jD,iD) + (1-omega)*T_prev(jD,iD);    
            % convergence criterion
        error = max(max(abs(T_diff - T_old)));
            % updating old values
        T_old = T_diff;
        n_iter = n_iter +1;
    end
        % updating previous values
        %\T_prev = T_diff;
        
      
        % Enforce BC outside the function
 
end