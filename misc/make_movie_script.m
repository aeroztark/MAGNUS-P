
% generic script to make movie from contour plots

loops = length(T_arr);
M(loops) = struct('cdata',[],'colormap',[]);

for j = 1:loops
    
    % enter the parameter here
contourf(X_KM,Z_KM,U(:,:,j).*SCALING_FACTOR,50,'Edgecolor','none');
colorbar
xlabel('x (km)')
ylabel('z (km)')
title(['scaled w (m/s) at time t=',num2str(T_arr(j)),'s'])
drawnow
M(j) = getframe(gcf);
end

% to play the movie later:
%movie(M,1,10)