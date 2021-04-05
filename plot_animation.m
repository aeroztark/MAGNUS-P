
% Script to launch animation of the simulation results

figure(1);
cmocean balance
Nmax=size(T_arr,2);


for n=1:Nmax
    
    %plot results: Pressure
    subplot(2,2,1);
    contourf(X_KM,Z_KM,P_PERT(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['Pressure Perturbation (Pa) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Temperature
    subplot(2,2,2);
    contourf(X_KM,Z_KM,T_PERT(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['Temperature Perturbation (K) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Horizontal Velocity
    subplot(2,2,3);
    contourf(X_KM,Z_KM,U(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['u (m/s) at time t=',num2str(T_arr(n)),'s']);
    
    %plot results: Vertical Velocity
    subplot(2,2,4);
    contourf(X_KM,Z_KM,W(:,:,n),50,'Edgecolor','none'); 
    colorbar; 
    xlabel('x (km)'); ylabel('z (km)'); title(['w (m/s) at time t=',num2str(T_arr(n)),'s']);
    
    pause(0.01)
end