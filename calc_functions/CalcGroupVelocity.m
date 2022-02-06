function [Cgx,Cgz] = CalcGroupVelocity(Lx,Lz,omega,u_wind,N)
% Function to calculate group velocity in x and z directions. Assumes
% Coriolis parameter ~ 0. inputs must be scalar.

% If only Cgz is needed, supply u_wind and N as 0

% Lx - horiz. wavelength in meters
% Lz - vertical wavelength in meters
% omega - intrinsic freq in rad/s
% u_wind - background u in m/s
% Cgx, Gz -> horizontal & vertical group velocities (m/s)


k = 2*pi/Lx;
m = 2*pi/Lz;
denom = omega*(k^2 + m^2);

Cgx = u_wind + k*(N^2-omega^2)/denom;

Cgz = -m*(omega^2)/denom;

end
