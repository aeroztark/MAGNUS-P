global wind
u_pert = U - wind(3:LastDomainZindex,3:end-2);
w_pert = W; % since no initial vertical velocity. If there is, need to subtract background vertical wind

f1 = u_pert.*w_pert; % this is 3d array of u'w'

f = sum(f1,2); % summed along x
f = f./abs(Xmax-Xmin);  % averaged by domain x width ( = Lx if periodic domain)
f = squeeze(f); % to make it 2d array (height vs time)

% contour plot
figure
contourf(T_arr,Z_KM,f,50,'Edgecolor','none')
xlabel('time (s)')
ylabel('z (km)')
title('Momentum Flux (m/s^2)')
colorbar

% vertical profile at a given time
% figure
% plot(f(:,165),Z_KM)
% xlabel('Momentum Flux')
% ylabel('z (km)')