Tprofile = 0;   % Enter Mars or Venus vertical T profile on z = 0:0.5:200 km
Lz = 1:0.5:50;  % km (Lz must be within 1-500 km for coeffs validity)
k = 100; % horizontal wavelength (km)

% supply the radiative coefficients file
coeffs = readtable('CO2dampingcoeffs.xlsx');

% -----------
m = (2*pi)./Lz; %(km^-1)

psi = log(m) - log(2*pi/500);

[nr,~] = size(coeffs);
tau_1ref = zeros(nr,length(psi)); %(tau_1 is damping rate while tau is timescale)

for i = 1:nr
    tau_1ref(i,:) = exp(coeffs.a(i) + psi.*coeffs.b(i) + (psi.^2).*coeffs.c(i) + (psi.^3).*coeffs.d(i)); % rates
end

% plotting for reference temperature profile:
contourf(Lz,coeffs.zref,1./tau_1ref)
xlabel('Vertical wavelength (km)')
ylabel('height (km)')
title('Radiative damping timescale (days)')
colorbar


% tau_1ref corresponds to Tref values in coefficients file. For application
% to any other temperature profile, interpolation is needed.The other
% temperature profile needs to be first interpolated to z array 0:0.5:200
% km, then scaling is applied
%tau_1 = tau_1ref.*(scaling(Tprofile)./scaling(coeffs.Tref));




function[dtheta_dt] = scaling(T)
% This function returns scaling parameter (dtheta/dt for any temperature
% array on vertical grid 0:0.5:200 km

dtheta_dt = (1.23686e5).*(exp(971./T)./((T.^2).*((exp(971./T)-1).^2)));

end