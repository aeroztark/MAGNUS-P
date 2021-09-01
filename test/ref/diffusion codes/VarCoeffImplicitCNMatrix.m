
function[A] = VarCoeffImplicitCNMatrix(k,dxz,i,j)
% This function generates the matrix for solution of 2D diffusion equation
% by Direct implicit method using Crank-Nicholson scheme.
% This function works for uniform grid (dx=dz) and VARIABLE diffusion
% coefficient, which applies to atmospheric physics problems.

% k -> diffusion constant vector of size jx1 (one k value for every j)
% dxz -> grid size for uniform grid (dx = dz) in meters
% i -> number of x values in u(length of x dimension)
% j -> number of z values in u (length of z dimension)

ij = i*j;

% create template arrays with ones and zeros
main_diag_temp = ones(ij,1);
plusi_diag_temp = ones(ij,1);
mini_diag_temp = plusi_diag_temp;

c_temp = ones(i,1);
min1c_unit = c_temp;
min1c_unit(end) = 0;
min1_diag_temp = repmat(min1c_unit,j,1);

plus1c_unit = c_temp;
plus1c_unit(1) = 0;
plus1_diag_temp = repmat(plus1c_unit,j,1);

% putting k into an (ij x 1) size vector (i.e. repeating each k value i
% times)
k_ij = repelem(k,i);
% make a and c arrays
a = 1+(2*k_ij./(dxz^2));
c = -k_ij./(2*dxz^2);

% make the diagonals
main_diag = main_diag_temp.*a;
plusi_diag = plusi_diag_temp.*c;
mini_diag = mini_diag_temp.*c;
plus1_diag = plus1_diag_temp.*c;
min1_diag = min1_diag_temp.*c;

A = spdiags([mini_diag,min1_diag,main_diag,plus1_diag,plusi_diag],[-i,-1,0,1,i],ij,ij);

end

