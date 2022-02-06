
% generic script to make movie from contour plots

loops = length(T_arr);
M(loops) = struct('cdata',[],'colormap',[]);

figure

for j = 1:loops
    
    % enter the parameter here

contourf(X_KM,Z_KM,W(:,:,j).*SCALING_FACTOR,50,'Edgecolor','none');
colorbar
xlabel('x (km)')
ylabel('z (km)')
title(['scaled w (m/s) at time t=',num2str(T_arr(j)),'s'])
drawnow
M(j) = getframe(gca);
end

% to playback:
% movie(gcf,M,1,0.5)

% % Save video
% v = VideoWriter('D:\address\filename.avi');
% open(v)
% writeVideo(v,M)
% close(v)