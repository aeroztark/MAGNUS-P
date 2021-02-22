function[T,rho,p,R,gamma,kvisc,H,C] = Earth_isothermal(z_array)

z_array = z_array';  % so that everything turns out as column vecotrs

T0 = 239;
p0 = 82311.6;
rho0 = 1.2;
H=p0/(rho0*9.8);
gamma = 1.4;
R = 287;
C = sqrt(gamma*R*T0);
N = sqrt(gamma-1)*9.8/C;

p = p0*exp(-z_array./H);
rho = rho0*exp(-z_array./H);
T = T0.*ones(length(z_array),1);
R = R.*ones(length(z_array),1);
C = C.*ones(length(z_array),1);
H = H.*ones(length(z_array),1);
gamma = gamma.*ones(length(z_array),1);

kvisc = (3.5e-7)*(T0^0.69)./rho; % Banks & Kockarts, 1973