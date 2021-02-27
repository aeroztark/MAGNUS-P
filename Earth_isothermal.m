function[T,rho,p,R,gamma,kvisc,H,C] = Earth_isothermal(z_array)

% Earth isothermal model to provide atmospheric state
% input: z_array -> Z x X array 
% all outputs have the same dimension as z_array

T_surface = 239;
p_surface = 1E5;
rho_surface = 1.2;
gamma = 1.4;    %constant gamma
R = 287;

H = p_surface/(rho_surface*9.8); % scaled height

p = p_surface*exp(-z_array./H);
rho = rho_surface*exp(-z_array./H);
T = T_surface.*ones(size(z_array));
R = R.*ones(size(z_array));
C = sqrt(gamma*R.*T);
C = C.*ones(size(z_array));
H = H.*ones(size(z_array));
gamma = gamma.*ones(size(z_array));
%N = sqrt(gamma-1)*9.8/C;
kvisc = (3.5e-7)*(T.^0.69)./rho; % Banks & Kockarts, 1973