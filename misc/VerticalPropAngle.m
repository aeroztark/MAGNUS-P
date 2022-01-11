function[beta] = VerticalPropAngle(k,m)

% Returns angle from vertical for high frequency gravity wave based on
% wavenumber geometry

% k, m are wavenumbers. Either both vectors of same length or one vector
% and the other a scalar

    beta = acosd(k./sqrt(k.^2 + m.^2));
    
end

%% Example script:
% Lx = linspace(0,500,500);   % 0-500 km, 500 pts
% Lz = linspace(0,50,500);    %0-50 km, 500 pts
% k = (2*pi)./Lx;
% m = (2*pi)./Lz;
% [K,M] = meshgrid(k,m);
% figure
% BETA = VerticalPropAngle(K,M);
% contourf(BETA);
% contourf((2*pi)./K,(2*pi)./M,BETA)
% xlabel('Lx')
% ylabel('Lz')
% title('Angle from vertical (deg)')
% colorbar