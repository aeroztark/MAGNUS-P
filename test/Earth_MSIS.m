
function[T,rho,p,R,gamma,kvisc,thermdiffus,H,C,omega_BV,omega_AC] = Earth_MSIS(z_array,lat,lon,year,DOY,UTsec)
% This function used NRLMSISE-00 model to compute background T, rho, P, gamma, R, kvisc, H, sound speed and BV/acoustic cutoff frequencies for Earth's atmosphere 
% input: provide z in meters and other params for MSISE (lat-lon & time)

warning('off')
[T,den] = atmosnrlmsise00(z_array(:,1),lat,lon,year,DOY,UTsec);

T = T(:,2); %ignore first row -> exospheric temperature
rho = den(:,6);  % 6th column is total mass density in kg/m3

% number densities of major species
N2 = den(:,3);
O2 = den(:,4);
O = den(:,2);

M = N2 + O2 + O;    % total number density

MM = (28.*N2 + 32.*O2 + 16.*O)./M; % mean molecular mass

Pr = 0.7;

% --- property computations
gamma = ((1.4.*(O2+N2))+1.67.*O)./M; % from Snively and Pasko. Simple weighing

kvisc = (((3.43.*N2 + 4.03.*O2 + 3.9.*O)./M)./rho).*(T.^0.69).*1e-7; % Rees 1989 % Eddy diffusion can be parameterized below 110 km

R = 8314./MM; % varying R with composition (in standard J/K/kg units), else constant value 287 can be used

p = rho.*R.*T; % using ideal gas equation to calculate pressure

H = R.*T./9.8;  % scale height

C = sqrt(gamma.*R.*T);  % speed of sound

omega_BV = sqrt(gamma-1).*9.8./C; % Brunt-Vaisala frequency

omega_AC = gamma.*9.8./(2.*C); % Acoustic cutoff frequency

% stretching column vectors to the same size as Z
[~,B] = size(z_array);
T = repmat(T,1,B); 
rho = repmat(rho,1,B);
p = repmat(p,1,B); 
R = repmat(R,1,B);
gamma = repmat(gamma,1,B); 
kvisc = repmat(kvisc,1,B); 
thermdiffus = kvisc./Pr;
H = repmat(H,1,B); 
C = repmat(C,1,B); 
omega_BV = repmat(omega_BV,1,B); 
omega_AC = repmat(omega_AC,1,B); 