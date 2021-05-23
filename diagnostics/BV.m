
function BV = BV(g,c_s,gamma)

% Compute Brunt-Vaisala frequency
% Inputs: 1D vertical arrays of g, speed of sound and gamma

%c_s = sqrt(gamma*R*T); % if calculating speed of sound here
g = [g(1,1); g];    % match other dimensions
BV = sqrt(gamma-1).*(g./c_s);   %in rad/s (divide by 2pi to get in Hz)

% figure
% plot(BV,z_c./1000)
% xlabel('BV (rad/s)')
% ylabel('h (km)')

end