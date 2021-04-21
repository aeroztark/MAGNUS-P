function[T_sat] = CO2_saturation(p)

% intpue -> p array (vertical pressure) and output -> T_sat array

% Saturation temperature for CO2 cloud formation on Mars

T0 = 136.3;
p0 = 100;
Lsv = 5.9e5;
R = 191.8;

T_sat = 1./((1/T0) - (R/Lsv).*log((1.35.*p)./p0));
end

function[] = plotCO2_clouds(T_sat,T)
% input -> T_sat -> saturation T array; T -> GW modulated T as Z by X
[row, col] = size(T);
Prob_cloud = zeros(size(T));

for i = 1:col
    Prob_cloud(:,i) = (T(:,i) < T_sat);
end

figure
contourf(Prob_cloud)

end

