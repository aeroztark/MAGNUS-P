function[w] = GaussianHillForcing(t,x_c,dx)

global domain

% x_c in meters
% dx in meters
H = 500; % meters (For hill, H>0 and for trough H<0)
L = 7000; % meters
u = 10; % m/s (can be +ve or -ve)
ramprate = (10/3600); % 0 to 10 m/s in 1 hr

topo = H.*exp(-(domain.x_c./L).^2); % meters
dh_dx = diff(topo)./domain.dx; % gradient
dh_dx = [dh_dx(1), dh_dx]; % append the same at position 1 to make the same length vector
w = (min(1,ramprate*t)*u).*dh_dx; % linearly ramp up and then hold constant

end