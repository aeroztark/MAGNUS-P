function[Ri] = Calc_Ri(BV,LastDomainZindex,U_PERT_profile,dz_m)

% compute Richardson Number for indication of wave breaking (Ri<0.25 ->
% shear instability & Ri<0 -> convective instability)

% Inputs:
% BV -> 1D vertical array of Brunt-Vaisala freq (rad/s)
% U_PERT_profile -> 1D vertical array of u perturbations (m/s)
% dz_m -> vertical grid spacing in meters

% Output -> 1D vertical array of Ri

den = (gradient(U_PERT_profile,dz_m)).^2;

Ri = (BV(3:LastDomainZindex).^2)./den;