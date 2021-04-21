% Not a working function. Just an FYI only

function[w] = MarsTopography_forcing(x_c)
% Use grabit.m to get points out as data
% [raw_x,raw_h]
% perform Savitzy-Golay fit to smoothen the profiles
%smooth_h = smoothdata(raw_h,'sgolay',3);

load('Smooth_pass260_topo.mat')

% sample at model x scale
hq = interp1(x,h,x_c./1000,'spline');

% Forcing from the obtained topography

dh_dx = diff(hq)./(dx/1000); % gradient
dh_dx = [dh_dx(1), dh_dx]; % append the same at position 1 to make the same length vector

% calculate forcing
u0 = 10;
w = u0.*dh_dx;  % surface wind x topography horizontal gradient

end